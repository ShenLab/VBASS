import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from torch.nn import init
from itertools import repeat
from torch.distributions import Normal


class GaussianSampler(nn.Module):
    """
    Layer that represents a sample from a
    Gaussian distribution.

    :param in_features: dimensions of the hidden layer.
           out_features: dimensions of the latent dim.
    """
    def __init__(self, in_features, out_features):
        super(GaussianSampler, self).__init__()
        self.in_features = in_features
        self.out_features = out_features

        self.mu = nn.Linear(in_features, out_features)
        self.sigma = nn.Linear(in_features, out_features)

    def forward(self, x):
        mu = self.mu(x)
        sigma = F.softplus(self.sigma(x))

        return Normal(mu, sigma).rsample(), mu, sigma


class SharedLayers(nn.Module):
    """
    Inference network shared part of Encoder and Classifier

    :param dims: dimensions of the networks
       given by the number of neurons on the form
       [input_dim, [hidden_dims], latent_dim].
    """
    def __init__(self,
                 dims,
                 use_bias: bool = True,
                 use_batch_norm: bool = False,
                 use_layer_norm: bool = True,
                 use_activation: bool = True,
                 activation_fn: nn.Module = nn.ELU,
                 dropout_rate: float = 0.1,
                 ):
        super(SharedLayers, self).__init__()

        [x_dim, h_dim] = dims
        neurons = [x_dim, *h_dim]
        linear_layers = []
        for i, (n_in, n_out) in enumerate(zip(neurons[:-1], neurons[1:])):
            linear_layers += [nn.Linear(in_features=n_in, out_features=n_out, bias=use_bias),
                              nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001) if use_batch_norm else None,
                              nn.LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
                              activation_fn() if use_activation else None,
                              nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None]
        linear_layers = [x for x in linear_layers if x]
        self.hidden = nn.Sequential(*linear_layers)

    def forward(self, x):
        return self.hidden(x)


class Encoder(nn.Module):
    def __init__(self,
                 dims,
                 sample_layer=GaussianSampler,
                 use_bias: bool = True,
                 use_batch_norm: bool = False,
                 use_layer_norm: bool = True,
                 use_activation: bool = True,
                 activation_fn: nn.Module = nn.ELU,
                 dropout_rate: float = 0.1,
                 ):
        """
        Inference network

        Attempts to infer the probability distribution
        p(z|x) from the data by fitting a variational
        distribution q_φ(z|x). Returns the two parameters
        of the distribution (µ, log σ²).

        :param dims: dimensions of the networks
           given by the number of neurons on the form
           [input_dim, [hidden_dims], latent_dim].
        """
        super(Encoder, self).__init__()

        [x_dim, h_dim, z_dim] = dims
        if len(h_dim) != 0:
            neurons = [x_dim, *h_dim]
            linear_layers = []
            for i, (n_in, n_out) in enumerate(zip(neurons[:-1], neurons[1:])):
                linear_layers += [nn.Linear(in_features=n_in, out_features=n_out, bias=use_bias),
                                  nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001) if use_batch_norm else None,
                                  nn.LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
                                  activation_fn() if use_activation else None,
                                  nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None]
            linear_layers = [x for x in linear_layers if x]

            self.hidden = nn.Sequential(*linear_layers)
            self.sample = sample_layer(h_dim[-1], z_dim)
        else:
            self.hidden = nn.Sequential(*[])
            self.sample = sample_layer(x_dim, z_dim)

    def forward(self, x):
        return self.sample(self.hidden(x))


class Classifier(nn.Module):
    def __init__(self,
                 dims,
                 use_bias: bool = True,
                 use_batch_norm: bool = False,
                 use_layer_norm: bool = True,
                 use_activation: bool = True,
                 activation_fn: nn.Module = nn.ELU,
                 dropout_rate: float = 0.1,):
        """
        Inference network for classification

        Attempts to infer the probability distribution
        p(y|x) from the data by fitting a variational
        distribution q_φ(y|x). Returns the parameter π
        of the categorical distribution p(y|π).

        :param dims: dimensions of the networks
           given by the number of neurons on the form
           [input_dim, [hidden_dims], latent_dim].
        """
        super(Classifier, self).__init__()
        [x_dim, h_dim, y_dim] = dims
        if len(h_dim) != 0:
            neurons = [x_dim, *h_dim]
            linear_layers = []
            for i, (n_in, n_out) in enumerate(zip(neurons[:-1], neurons[1:])):
                linear_layers += [nn.Linear(in_features=n_in, out_features=n_out, bias=use_bias),
                                  nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001) if use_batch_norm else None,
                                  nn.LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
                                  activation_fn() if use_activation else None,
                                  nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None]
            linear_layers = [x for x in linear_layers if x]

            self.hidden = nn.Sequential(*linear_layers)
            self.logits = nn.Linear(h_dim[-1], y_dim)
        else:
            self.hidden = nn.Sequential(*[])
            self.logits = nn.Linear(x_dim, y_dim)

    def forward(self, x):
        pi = F.log_softmax(self.logits(self.hidden(x)), dim=-1)
        return pi


class Decoder(nn.Module):
    def __init__(self,
                 dims,
                 output_activation,
                 use_bias: bool = True,
                 use_batch_norm: bool = False,
                 use_layer_norm: bool = True,
                 use_activation: bool = True,
                 activation_fn: nn.Module = nn.ELU,
                 dropout_rate: float = 0,
                 ):
        """
        Generative network

        Generates samples from the original distribution
        p(x) by transforming a latent representation, e.g.
        by finding p_θ(x|z).

        :param dims: dimensions of the networks
            given by the number of neurons on the form
            [latent_dim, [hidden_dims], input_dim].
        """
        super(Decoder, self).__init__()

        [z_dim, h_dim, x_dim] = dims

        neurons = [z_dim, *h_dim]
        linear_layers = []
        for i, (n_in, n_out) in enumerate(zip(neurons[:-1], neurons[1:])):
            linear_layers += [nn.Linear(in_features=n_in, out_features=n_out, bias=use_bias),
                              nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001) if use_batch_norm else None,
                              nn.LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
                              activation_fn() if use_activation else None,
                              nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None]
        linear_layers = [x for x in linear_layers if x]
        self.hidden = nn.Sequential(*linear_layers)

        self.reconstruction = nn.Linear(h_dim[-1], x_dim)

        self.output_activation = output_activation

    def forward(self, x, y=None):
        # suppose y is a boolean vector of whether return true value
        x = self.hidden(x)
        if y is not None:
            # needs poisson activation
            return self.output_activation(self.reconstruction(x), y)
        else:
            return self.output_activation(self.reconstruction(x))


class Decoder_v2(nn.Module):
    def __init__(self,
                 dims,
                 out_likelihood,
                 is_null_hypothesis,
                 use_bias: bool = True,
                 use_batch_norm: bool = False,
                 use_layer_norm: bool = True,
                 use_activation: bool = True,
                 activation_fn: nn.Module = nn.ELU,
                 dropout_rate: float = 0,
                 ):
        """
        Generative network

        Generates samples from the original distribution
        p(x) by transforming a latent representation, e.g.
        by finding p_θ(x|z).

        :param dims: dimensions of the networks
            given by the number of neurons on the form
            [latent_dim, [hidden_dims], input_dim].
        """
        super(Decoder_v2, self).__init__()

        [z_dim, h_dim, x_dim] = dims

        self.is_null_hypothesis = is_null_hypothesis
        self.out_likelihood = out_likelihood
        self.x_dim = x_dim

        if self.is_null_hypothesis and (self.out_likelihood in ['Poisson', 'NegativeBinomial']):
            self.hidden = None
            self.reconstruction = None
        else:
            neurons = [z_dim, *h_dim]
            linear_layers = []
            for i, (n_in, n_out) in enumerate(zip(neurons[:-1], neurons[1:])):
                linear_layers += [nn.Linear(in_features=n_in, out_features=n_out, bias=use_bias),
                                  nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001) if use_batch_norm else None,
                                  nn.LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
                                  activation_fn() if use_activation else None,
                                  nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None]
            linear_layers = [x for x in linear_layers if x]
            self.hidden = nn.Sequential(*linear_layers)
            # for binomial distribution, output beta distribution parameters
            # for poisson distribution, output gamma distribution parameters
            # so, anyway, there will be two reconstructions, but I suppose activation will be same
            self.reconstruction_1 = nn.Sequential(*[nn.Linear(h_dim[-1], x_dim),
                                                    nn.Softplus()])
            self.reconstruction_2 = nn.Sequential(*[nn.Linear(h_dim[-1], x_dim),
                                                    nn.Softplus()])

    def forward(self, x):
        # suppose y is a boolean vector of whether return true value
        if self.hidden is None:
            return torch.zeros_like(x, device=x.device)
        else:
            x = self.hidden(x)
            return self.reconstruction_1(x), self.reconstruction_2(x)


class PassThroughActivation(nn.Module):
    def __init__(self):
        super(PassThroughActivation, self).__init__()

    def forward(self, x):
        return x


class PoissonActivation(nn.Module):
    def __init__(self):
        super(PoissonActivation, self).__init__()
        self.activation = nn.Softplus()

    def forward(self, x):
        result = self.activation(x)
        # result += torch.ones_like(result)
        # result[y, :] = result[y, :].detach()
        # result.register_hook(lambda grad: grad * y.float())
        return result



