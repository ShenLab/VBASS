from typing import List
from typing_extensions import Literal

from .module import *
from .utils import *
from torch.distributions import Normal, Poisson, Categorical
from torch.distributions import kl_divergence


class DeepGenerativeMixModel(nn.Module):
    def __init__(self,
                 x_in_dim: List[int],
                 h_in_shared_dim: List[int] = None,
                 h_in_z_dim: List[int] = None,
                 h_in_y_dim: List[int] = None,
                 z_dim: int = None,
                 y_dim: int = None,
                 y_kl_param: List[float] = None,
                 y_kl_param_learnable: bool = False,
                 gamma_kl_param_learnable: bool = True,
                 gamma_alpha_kl_param: List[float] = None,
                 gamma_beta_kl_param: List[float] = None,
                 gamma_beta_based_on_alpha: bool = False,
                 gamma_fixed: bool = False,
                 h_out_shared_dim: List[int] = None,
                 h_out_dims: List[List[int]] = None,
                 x_out_dim: List[int] = None,
                 x_out_likelihood: List[Literal['Binomial', 'Gaussian', 'Poisson', 'NegativeBinomial']] = None,
                 x_out_mixture: List[bool] = None,
                 null_hypothesis: int = 0,
                 batch_norm: List[bool] = None,
                 layer_norm: List[bool] = None,
                 ):
        """
        mixture model of TADA inspired from the paper
        'Semi-Supervised Learning with Deep Generative Models'
        (Kingma 2014) in PyTorch.

        This "Semi-supervised model" is a probabilistic
        model that incorporates label information in both
        inference and generation.

        Initialise a new generative model
        :param x_in_dim: array-like, dims of different input
        :param h_in_shared_dim: array-like, shared encoder
        :param h_in_z_dim: array-like, z-specific encoder, ended with a sample layer
        :param h_in_y_dim: array-like, y-specific encoder, ended with softmax layer
        :param z_dim: array-like, z dimensions
        :param y_dim: array-like, one-hot classifier
        :param h_out_shared_dim: shared decoder dims.
        :param h_out_dims: decoder dims.
        :param x_out_likelihood: array-like, output likelihood
        :param x_out_mixture: array-like, boolean, whether to use mixture out
        """
        # set default parameters
        h_in_shared_dim = [32] if h_in_shared_dim is None else list(h_in_shared_dim)
        h_in_z_dim = [] if h_in_z_dim is None else list(h_in_z_dim)
        h_in_y_dim = [] if h_in_y_dim is None else list(h_in_y_dim)
        z_dim = 8 if z_dim is None else int(z_dim)
        y_dim = 2 if y_dim is None else int(y_dim)
        y_kl_param = [1.25, -1.25] if y_kl_param is None else list(y_kl_param)
        gamma_alpha_kl_param = [16.61025718, 6.77746547] if gamma_alpha_kl_param is None else list(gamma_alpha_kl_param)
        gamma_beta_kl_param = [1.0, 1.0] if gamma_beta_kl_param is None else list(gamma_beta_kl_param)
        h_out_shared_dim = [] if h_out_shared_dim is None else list(h_out_shared_dim)
        x_out_dim = x_in_dim.copy() if x_out_dim is None else list(x_out_dim)
        h_out_dims = [[32]] * len(x_out_dim) if h_out_dims is None else list(h_out_dims)
        x_out_mixture = [False, True] if x_out_mixture is None else list(x_out_mixture)
        batch_norm = [False, False] if batch_norm is None else list(batch_norm)
        layer_norm = [False, False] if layer_norm is None else list(layer_norm)

        if x_out_likelihood is None:
            raise ValueError('likelihood should be specified')

        # set out decoder dims

        self.x_in_dim = x_in_dim
        self.y_dim = y_dim
        self.x_out_dim = x_out_dim
        self.x_out_likelihood = x_out_likelihood
        self.x_out_mixture = x_out_mixture
        self.x_mixture_parts = [i for i in range(len(self.x_out_mixture)) if self.x_out_mixture[i]]
        self.x_unmixture_parts = [i for i in range(len(self.x_out_mixture)) if not self.x_out_mixture[i]]

        super(DeepGenerativeMixModel, self).__init__()

        # Encoder has two parts, first a shared part
        self.shared_encoder = SharedLayers([sum(x_in_dim), h_in_shared_dim],
                                           use_batch_norm=batch_norm[0],
                                           use_layer_norm=layer_norm[0])
        encoder_shared_dim = h_in_shared_dim[-1] if len(h_in_shared_dim) != 0 else sum(x_in_dim)
        # Then unique parts for z and y
        self.encoder = Encoder([encoder_shared_dim, h_in_z_dim, z_dim],
                               use_batch_norm=batch_norm[0],
                               use_layer_norm=layer_norm[0])
        # in this model, classifier is placed after z.
        self.classifier = Classifier([encoder_shared_dim, h_in_y_dim, self.y_dim],
                                     use_batch_norm=batch_norm[0],
                                     use_layer_norm=layer_norm[0])
        # Only one decoder, takes both z and y as input, note x_out does not have to be original x
        self.shared_decoder = SharedLayers([z_dim, h_out_shared_dim],
                                           use_batch_norm=batch_norm[1],
                                           use_layer_norm=layer_norm[1])
        decoder_shared_dim = h_out_shared_dim[-1] if len(h_out_shared_dim) != 0 else z_dim
        self.null_hypothesis = null_hypothesis

        decoder_list = []
        for i in range(len(self.x_out_mixture)):
            if self.x_out_mixture[i]:
                # if mixture output
                mixture_decoders = []
                for k in range(self.y_dim):
                    # for each hypothesis
                    mixture_decoders.append(
                        Decoder_v2([decoder_shared_dim,
                                    h_out_dims[i],
                                    x_out_dim[i]],
                                   x_out_likelihood[i],
                                   k == self.null_hypothesis,
                                   use_batch_norm=batch_norm[1],
                                   use_layer_norm=layer_norm[1])
                    )
                decoder_list.append(nn.ModuleList(mixture_decoders))
            else:
                # only need one decoder
                decoder_list.append(
                    Decoder_v2([decoder_shared_dim,
                                h_out_dims[i],
                                x_out_dim[i]],
                               x_out_likelihood[i],
                               False,
                               use_batch_norm=batch_norm[1],
                               use_layer_norm=layer_norm[1])
                )
        self.decoders = nn.ModuleList(decoder_list)

        # set global-wise y categorical distribution parameters
        self.y_kl_param = torch.nn.Parameter(torch.zeros(y_dim),
                                             requires_grad=y_kl_param_learnable) \
            if y_kl_param is None else torch.nn.Parameter(torch.tensor(y_kl_param),
                                                          requires_grad=y_kl_param_learnable)

        self.gamma_fixed = gamma_fixed
        self.gamma_beta_based_on_alpha = gamma_beta_based_on_alpha

        self.gamma_alpha_kl_param = torch.nn.Parameter(torch.tensor(gamma_alpha_kl_param),
                                                       requires_grad=gamma_kl_param_learnable)
        if not self.gamma_beta_based_on_alpha:
            self.gamma_beta_kl_param = torch.nn.Parameter(torch.tensor(gamma_beta_kl_param),
                                                          requires_grad=gamma_kl_param_learnable)
        else:
            self.gamma_beta_kl_param = None
            

        self.kl_z = None
        self.kl_y = None

    def inference(self, x):
        # generate latent variable
        encoder_shared = self.shared_encoder(x)
        # generate z and reparameterize
        z, z_mu, z_sigma = self.encoder(encoder_shared)
        # y is generated after z
        y_logits = self.classifier(encoder_shared)

        self.kl_z = self._kld_z((z_mu, z_sigma))
        self.kl_y = self._kld_y(y_logits,
                                self.y_kl_param.repeat(y_logits.shape[0], 1))
        return z, z_mu, z_sigma, y_logits

    def _get_gamma_beta_kl_param(self):
        if self.gamma_beta_based_on_alpha:
            return torch.exp(6.7771073*self.gamma_alpha_kl_param**-1.7950864 - 0.2168248)
        else:
            return self.gamma_beta_kl_param

    def generative(self, z, y):
        # assume y is only single value that indicate the hypothesis.
        decoder_shared = self.shared_decoder(z)
        x_mu = []
        for i in range(len(self.decoders)):
            if self.x_out_likelihood[i] == 'Poisson' and self.gamma_fixed:
                x_mu.append((self.gamma_alpha_kl_param.repeat(decoder_shared.shape[0], 1),
                             self._get_gamma_beta_kl_param().repeat(decoder_shared.shape[0], 1)))
            else:
                if self.x_out_mixture[i]:
                    x_mu.append(self.decoders[i][y](decoder_shared))
                else:
                    x_mu.append(self.decoders[i](decoder_shared))
        return x_mu

    def forward(self, x, y):
        # assume y is only single value that indicate the hypothesis.
        z, z_mu, z_sigma, y_logits = self.inference(x)
        # Reconstruct data point from latent data and label
        x_mu = self.generative(z, y)

        return x_mu, z_mu, z_sigma, y_logits

    def _kld_z(self, q_param, p_param=None):
        """
        Computes the KL-divergence of
        some element z.

        KL(q||p) = -∫ q(z) log [ p(z) / q(z) ]
                 = -E[log p(z) - log q(z)]

        :param z: sample from q-distribution
        :param q_param: (mu, log_var) of the q-distribution
        :param p_param: (mu, log_var) of the p-distribution
        :return: KL(q||p)
        """
        (mu, sigma) = q_param
        qz = Normal(mu, sigma)
        if p_param is None:
            pz = Normal(torch.zeros_like(mu), torch.ones_like(sigma))
        else:
            (mu, sigma) = p_param
            pz = Normal(mu, sigma)
        result = kl_divergence(qz, pz).sum(-1)
        return result

    def _kld_y(self, q_param, p_param=None):
        """
        Computes the KL-divergence of
        some element y.

        KL(q||p) = -∫ q(y) log [ p(y) / q(y) ]
                 = -E[log p(y) - log q(y)]

        :param z: sample from q-distribution
        :param q_param: pi of the q-distribution Categorical
        :param p_param: pi of the p-distribution Categorical
        :return: KL(q||p)
        """
        qy = Categorical(logits=q_param)

        if p_param is None:
            py = Categorical(logits=torch.ones_like(q_param))
        else:
            py = Categorical(logits=p_param)
        result = kl_divergence(qy, py)
        return result

    def sample(self, z_mu, z_sigma, y_logits):
        """
        Samples from the Decoder to generate an x.
        :param z_mu: latent normal variable mean
        :param z_sigma: latent normal variable std
        :param y_logits: label (one-hot encoded)
        :return: x
        """
        z = Normal(z_mu, z_sigma).rsample()
        y = Categorical(logits=y_logits).rsample()
        x = self.generative(z, y)
        return x
