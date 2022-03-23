import numpy as np
import torch.distributions
from pyro.distributions.conjugate import BetaBinomial
from pyro.distributions.conjugate import GammaPoisson
from torch.distributions import kl_divergence
from .model import *


class LikelihoodMixModel(nn.Module):
    def __init__(self, model: DeepGenerativeMixModel):
        """
        A helper Class to calculate likelihood loss
        Support Poisson and Binomial distributions for now
        TODO:
        Add Negative binomial and Gaussian Distributions
        :param model: A DeepGenerativeMixModel_v2
        """
        super(LikelihoodMixModel, self).__init__()
        self.x_split = model.x_out_dim
        self.x_out_likelihood = model.x_out_likelihood
        self.null_hypothesis = model.null_hypothesis
        self.model = model

    def forward(self, xs, reconstruction, x_vars, y):
        xs_split = torch.split(xs, self.x_split, dim=1)
        x_vars_split = torch.split(x_vars, self.x_split, dim=1)
        likelihood = []
        penalty = []
        for i in range(len(self.x_out_likelihood)):
            if self.x_out_likelihood[i] == "Binomial":
                likelihood.append(
                    BetaBinomial(
                        concentration1=reconstruction[i][0],
                        concentration0=reconstruction[i][1],
                        total_count=x_vars_split[i]
                    ).log_prob(xs_split[i]).sum(dim=1, keepdims=True)
                )
            elif self.x_out_likelihood[i] == "Poisson":
                if (y == self.null_hypothesis):
                    likelihood.append(
                        torch.distributions.Poisson(
                            rate=x_vars_split[i]*2
                        ).log_prob(xs_split[i]).sum(dim=1, keepdims=True)
                    )
                else:
                    likelihood.append(
                        GammaPoisson(
                            concentration=reconstruction[i][0]*x_vars_split[i]*2,
                            rate=reconstruction[i][1],
                        ).log_prob(xs_split[i]).sum(dim=1, keepdims=True)
                    )
                    penalty.append(
                        # if it is alter hypothesis, it cannot be too close to 0
                        kl_divergence(
                            torch.distributions.Gamma(
                                reconstruction[i][0],
                                reconstruction[i][1]
                            ),
                            torch.distributions.Gamma(
                                self.model.gamma_alpha_kl_param.repeat(reconstruction[i][0].shape[0], 1),
                                self.model.gamma_beta_kl_param.repeat(reconstruction[i][1].shape[0], 1)
                            )
                        ).sum(dim=1, keepdims=True)
                    )
        if (y == self.null_hypothesis):
            return torch.cat(likelihood, dim=1), None
        else:
            return torch.cat(likelihood, dim=1), torch.cat(penalty, dim=1).sum(dim=1)


def ELBO_mixmodel(model, x_in, x_out, x_var, y, likelihood_weight,
                  device, writer, counter,
                  alpha=repeat(1), beta=repeat(1), gamma=repeat(1), epsilon=repeat(1)):
    # alpha is classification loss
    # beta is kl loss
    # gamma is penalty loss
    likelihood_fn = LikelihoodMixModel(model)
    is_labelled = ~torch.isnan(y).resize(y.shape[0],)

    # Prepare for sampling
    x_in_s_labeled = x_in[is_labelled, ]
    x_out_s_labeled = x_out[is_labelled, ]
    x_vars_labeled = x_var[is_labelled, ]
    ys_labeled = y[is_labelled, ]

    x_in_s_unlabeled = x_in[~is_labelled, ]
    x_out_s_unlabeled = x_out[~is_labelled, ]
    x_vars_unlabeled = x_var[~is_labelled, ]
    ys_unlabeled = y[~is_labelled, ]
    input_unlabeled = is_labelled[~is_labelled]

    x_ins = x_in_s_labeled
    x_outs = x_out_s_labeled
    x_vars = x_vars_labeled
    ys = onehot(ys_labeled.long(), model.y_dim, device="cpu")
    labeled = ~torch.isnan(y)[is_labelled, ]

    for i in range(x_in_s_unlabeled.shape[0]):
        y_repeat = enumerate_discrete(x_in_s_unlabeled[i, ].resize(1, x_in_s_unlabeled.shape[1]), model.y_dim)
        x_in_repeat = x_in_s_unlabeled[i, ].repeat(model.y_dim, 1)
        x_out_repeat = x_out_s_unlabeled[i, ].repeat(model.y_dim, 1)
        x_var_repeat = x_vars_unlabeled[i, ].repeat(model.y_dim, 1)
        label_repeat = input_unlabeled[i].repeat(model.y_dim, 1)

        x_ins = torch.cat((x_ins, x_in_repeat))
        x_outs = torch.cat((x_outs, x_out_repeat))
        ys = torch.cat((ys, y_repeat))
        x_vars = torch.cat((x_vars, x_var_repeat))
        labeled = torch.cat((labeled, label_repeat))

    labeled = labeled.detach().cpu().numpy()[:, 0]
    x_ins = x_ins.to(device)
    x_outs = x_outs.to(device)
    x_vars = x_vars.to(device)
    ys = ys.to(device)
    labels = reverse_onehot(ys).detach().cpu().numpy()
    # calculate likelihood for each hypothesis separately
    likelihood_unmixture = []
    likelihood_mixture = []
    penalty = []
    kl_labeled_y = []
    kl_unlabeled_y = []
    kl_labeled_z = []
    kl_unlabeled_z = []
    classification_loss = []
    for i in range(model.y_dim):
        x_ins_i = x_ins[labels == i, :]
        x_outs_i = x_outs[labels == i, :]
        x_vars_i = x_vars[labels == i, :]
        labeled_i = labeled[labels == i]
        ys_i = ys[labels == i, :]
        reconstruction_i, _, _, y_logits_i = model(x_ins_i, i)
        # Add likelihood and penalty
        likelihood_i, penalty_i = likelihood_fn(x_outs_i, reconstruction_i, x_vars_i, i)
        likelihood_i = likelihood_i * likelihood_weight.to(device).repeat(likelihood_i.shape[0], 1)

        likelihood_i_mixture = likelihood_i[:, model.x_mixture_parts].sum(dim=1, keepdims=True)
        if len(model.x_unmixture_parts):
            likelihood_i_unmixture = likelihood_i[:, model.x_unmixture_parts].sum(dim=1, keepdims=True)
            likelihood_unmixture.append(likelihood_i_unmixture)

        likelihood_i_mixture_unlabeled = likelihood_i_mixture[~labeled_i, 0] + y_logits_i[ys_i.bool()][~labeled_i]
        likelihood_mixture.append(likelihood_i_mixture_unlabeled.resize(sum(~labeled_i), 1))
        likelihood_i_mixture_labeled = likelihood_i_mixture[labeled_i, ]
        likelihood_unmixture.append(likelihood_i_mixture_labeled)

        if penalty_i is not None:
            penalty_unlabeled_i = penalty_i[~labeled_i, ] * torch.exp(y_logits_i[~labeled_i, i])
            penalty_labeled_i = penalty_i[labeled_i, ]
            penalty.append(torch.cat((penalty_unlabeled_i, penalty_labeled_i)))
        # Add KL
        kl_i_y = model.kl_y
        kl_i_z = model.kl_z
        kl_labeled_y.append(kl_i_y[labeled_i, ])
        kl_unlabeled_y.append(kl_i_y[~labeled_i, ])
        kl_labeled_z.append(kl_i_z[labeled_i, ])
        kl_unlabeled_z.append(kl_i_z[~labeled_i, ])
        # Add classification loss
        y_logits_labeled_i = y_logits_i[labeled_i, ]
        ys_labeled_i = ys_i[labeled_i, ]
        if len(ys_labeled_i) > 0:
            classification_loss.append((y_logits_labeled_i * ys_labeled_i).sum(dim=1))
    likelihood_mixture = torch.logsumexp(torch.cat(likelihood_mixture, dim=1), dim=1, keepdim=True)
    likelihood_unmixture = torch.cat(likelihood_unmixture)
    likelihood = torch.cat((likelihood_mixture, likelihood_unmixture)).mean()
    kl_y = torch.cat((torch.cat(kl_labeled_y), kl_unlabeled_y[0])).mean()
    kl_z = torch.cat((torch.cat(kl_labeled_z), kl_unlabeled_z[0])).mean()
    classification_loss = torch.cat(classification_loss).mean() if len(classification_loss) else None
    penalty = torch.cat(penalty).mean()
    # ADD elbo
    elbo = -likelihood + next(beta) * kl_y + next(gamma) * kl_z + next(epsilon) * penalty
    # ADD writer
    writer.add_scalar('Likelihood/likelihood_mixture', likelihood_mixture.mean().item(), counter)
    writer.add_scalar('Likelihood/likelihood_unmixture', likelihood_unmixture.mean().item(), counter)
    writer.add_scalar('Likelihood/likelihood', likelihood.item(), counter)
    writer.add_scalar('Likelihood/elbo', elbo.item(), counter)
    writer.add_scalar('Likelihood/penalty', penalty.item(), counter)
    writer.add_scalar('Likelihood/kl_y', kl_y.item(), counter)
    writer.add_scalar('Likelihood/kl_z', kl_z.item(), counter)
    if classification_loss is not None:
        writer.add_scalar('Likelihood/classification', classification_loss.item(), counter)
    # ADD loss
    if classification_loss is not None:
        loss = elbo - next(alpha) * classification_loss
    else:
        loss = elbo
    return loss


def logBayesFactor_mixmodel(model, x_in, x_out, x_var, y, likelihood_weight,
                            device="cuda"):
    # x_in is the values to encoder
    # x_out_vars are the additional data
    likelihood_fn = LikelihoodMixModel(model)
    is_labelled = ~torch.isnan(y).resize(y.shape[0],)
    # Prepare for sampling
    x_in_s_labeled = x_in[is_labelled, ]
    x_out_s_labeled = x_out[is_labelled, ]
    x_vars_labeled = x_var[is_labelled, ]
    ys_labeled = y[is_labelled, ]
    x_in_s_unlabeled = x_in[~is_labelled, ]
    x_out_s_unlabeled = x_out[~is_labelled, ]
    x_vars_unlabeled = x_var[~is_labelled, ]
    ys_unlabeled = y[~is_labelled, ]
    input_unlabeled = is_labelled[~is_labelled]

    x_ins = x_in_s_labeled
    x_outs = x_out_s_labeled
    x_vars = x_vars_labeled
    ys = onehot(ys_labeled.long(), model.y_dim, device="cpu")
    labeled = ~torch.isnan(y)[is_labelled, ]

    for i in range(x_in_s_unlabeled.shape[0]):
        y_repeat = enumerate_discrete(x_in_s_unlabeled[i, ].resize(1, x_in_s_unlabeled.shape[1]), model.y_dim)
        x_in_repeat = x_in_s_unlabeled[i, ].repeat(model.y_dim, 1)
        x_out_repeat = x_out_s_unlabeled[i, ].repeat(model.y_dim, 1)
        x_var_repeat = x_vars_unlabeled[i, ].repeat(model.y_dim, 1)
        label_repeat = input_unlabeled[i].repeat(model.y_dim, 1)

        x_ins = torch.cat((x_ins, x_in_repeat))
        x_outs = torch.cat((x_outs, x_out_repeat))
        ys = torch.cat((ys, y_repeat))
        x_vars = torch.cat((x_vars, x_var_repeat))
        labeled = torch.cat((labeled, label_repeat))

    labeled = labeled.detach().cpu().numpy()[:, 0]
    x_ins = x_ins.to(device)
    x_outs = x_outs.to(device)
    x_vars = x_vars.to(device)
    ys = ys.to(device)
    labels = reverse_onehot(ys).detach().cpu().numpy()
    # calculate likelihood for each hypothesis separately
    likelihood_mixture = []
    reconstruction_mean = []
    reconstruction_var = []
    y_logits = None
    for i in range(model.y_dim):
        x_ins_i = x_ins[labels == i, :]
        x_outs_i = x_outs[labels == i, :]
        x_vars_i = x_vars[labels == i, :]
        labeled_i = labeled[labels == i]
        ys_i = ys[labels == i, :]
        # repeat 100 times
        x_ins_i = x_ins_i.repeat(100, 1)
        x_outs_i = x_outs_i.repeat(100, 1)
        ys_i = ys_i.repeat(100, 1)
        x_vars_i = x_vars_i.repeat(100, 1)
        labeled_i = labeled_i.repeat(100)
        reconstruction_i, z_mu_i, z_sigma_i, y_logits_i = model(x_ins_i, i)
        y_logits = torch.split(y_logits_i, int(sum(labels == i)))[0]
        z_mu = torch.split(z_mu_i, int(sum(labels == i)))[0]
        z_sigma = torch.split(z_sigma_i, int(sum(labels == i)))[0]

        # Add likelihood and penalty
        likelihood_i, _ = likelihood_fn(x_outs_i, reconstruction_i, x_vars_i, i)
        likelihood_i = likelihood_i * likelihood_weight.to(device).repeat(likelihood_i.shape[0], 1)
        likelihood_i_mixture = likelihood_i[:, model.x_mixture_parts]

        likelihood_i_mixture_unlabeled = likelihood_i_mixture[~labeled_i, 0] + y_logits_i[ys_i.bool()][~labeled_i]
        likelihood_mixture.append(
            torch.stack(torch.split(likelihood_i_mixture_unlabeled,
                                    int(sum(labels == i)))).mean(dim=0, keepdims=True).T
        )
        if i != model.null_hypothesis:
            reconstruction_mean.append(
                torch.stack(torch.split(reconstruction_i[1][0]/reconstruction_i[1][1],
                                        int(sum(labels == i)))).mean(dim=0)
            )
            reconstruction_var.append(
                torch.stack(torch.split(reconstruction_i[1][1],
                                        int(sum(labels == i)))).mean(dim=0)
            )
    likelihood_mixture = torch.cat(likelihood_mixture, dim=1)
    BF = likelihood_mixture[:, 1:] - likelihood_mixture[:, 0].repeat(likelihood_mixture.shape[1] - 1, 1).T

    return BF, y_logits, z_mu, z_sigma, torch.cat(reconstruction_mean), torch.cat(reconstruction_var)
