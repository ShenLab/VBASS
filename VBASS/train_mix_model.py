import sys
import os
import argparse
import json
from Model import *
import matplotlib.pyplot as plt
import time
from torch.utils.tensorboard import SummaryWriter


def run_one_epoch(counter, train_flag, dataloader, model, likelihood_weight,
                  optimizer, writer, debug_dir, device="cuda",
                  alpha=repeat(1), beta=repeat(1), gamma=repeat(1), epsilon=repeat(1)):
    torch.set_grad_enabled(train_flag)
    model.train() if train_flag else model.eval()

    losses = []

    for i, (x_in, x_out, x_var, x_label) in enumerate(dataloader):  # collection of tuples with iterator

        loss = ELBO_mixmodel(model, x_in, x_out, x_var, x_label,
                             likelihood_weight, device, writer, counter,
                             alpha=alpha,
                             beta=beta,
                             gamma=gamma,
                             epsilon=epsilon)
        if train_flag:
            if torch.isnan(loss):
                torch.save(x_var, f'{debug_dir}/x_var.nan.loss.{counter}.pt')
                torch.save(x_label, f'{debug_dir}/x_label.nan.loss.{counter}.pt')
                torch.save(model.state_dict(), f'{debug_dir}/model.nan.loss.{counter}.pt')
                break
            else:
                torch.save(x_var, f'{debug_dir}/x_var.{counter}.pt')
                torch.save(x_label, f'{debug_dir}/x_label.{counter}.pt')
                torch.save(model.state_dict(), f'{debug_dir}/model.{counter}.pt')

            loss.backward()  # back propagation

            writer.add_scalar('Loss/train', loss.item(), counter)
            writer.add_scalar('Parameter/y_kl[0]', model.y_kl_param[0].item(), counter)
            writer.add_scalar('Parameter/y_kl[1]', model.y_kl_param[1].item(), counter)
            writer.add_scalar('Parameter/gamma_alpha_kl_param[0]', model.gamma_alpha_kl_param[0].item(), counter)
            writer.add_scalar('Parameter/gamma_alpha_kl_param[1]', model.gamma_alpha_kl_param[1].item(), counter)
            writer.add_scalar('Parameter/gamma_beta_kl_param[0]', model.gamma_beta_kl_param[0].item(), counter)
            writer.add_scalar('Parameter/gamma_beta_kl_param[1]', model.gamma_beta_kl_param[1].item(), counter)
            for tag, value in model.named_parameters():
                tag = tag.replace('.', '/')
                writer.add_histogram('weights/'+tag, value.data.cpu().numpy(), counter)
                try:
                    writer.add_histogram('grads/'+tag, value.grad.data.cpu().numpy(), counter)
                except ValueError:
                    print(f"failed to add grad histogram for '{tag}' in counter: {counter}")
                except AttributeError:
                    if tag != "y_kl_param":
                        print(f"failed to add grad histogram for '{tag}' in counter: {counter}")

            optimizer.step()
            optimizer.zero_grad()
            counter += 1

            losses.append(loss.detach().cpu().numpy())
            print(f'finished batch {i} with {losses[-1]}')

    return np.mean(losses), counter


def output_log_BF(dataloader, model, likelihood_weight, device="cuda"):
    torch.set_grad_enabled(False)
    model.eval()

    logBFs = None
    reconstruction_means = None
    reconstruction_vars = None
    y_logits = None
    z_sigmas = None
    z_mus  =None

    for i, (x_in, x_out, x_var, x_label) in enumerate(dataloader):  # collection of tuples with iterator

        logBF, y_logit, z_mu, z_sigma, reconstruction_mean, reconstruction_var \
            = logBayesFactor_mixmodel(model, x_in, x_out, x_var, x_label, likelihood_weight, device)

        logBFs = np.concatenate((logBFs, logBF.detach().cpu().numpy())) \
            if logBFs is not None else logBF.detach().cpu().numpy()

        reconstruction_means = np.concatenate((reconstruction_means, reconstruction_mean.detach().cpu().numpy())) \
            if reconstruction_means is not None else reconstruction_mean.detach().cpu().numpy()

        reconstruction_vars = np.concatenate((reconstruction_vars, reconstruction_var.detach().cpu().numpy())) \
            if reconstruction_vars is not None else reconstruction_var.detach().cpu().numpy()

        y_logits = np.concatenate((y_logits, y_logit.detach().cpu().numpy())) \
            if y_logits is not None else y_logit.detach().cpu().numpy()

        z_mus = np.concatenate((z_mus, z_mu.detach().cpu().numpy())) \
            if z_mus is not None else z_mu.detach().cpu().numpy()

        z_sigmas = np.concatenate((z_sigmas, z_sigma.detach().cpu().numpy())) \
            if z_sigmas is not None else z_sigma.detach().cpu().numpy()

    return logBFs, y_logits, z_mus, z_sigmas, reconstruction_means, reconstruction_vars


def train_model(config):
    torch.set_num_threads(9)
    model_name = DeepGenerativeMixModel
    data_dir = config['input']['data_dir']
    model_dir = config['output']['model_dir']
    input_x = pd.read_csv(f'{data_dir}/input.x.csv', index_col=0)
    input_x_out = pd.read_csv(f'{data_dir}/input.x.csv', index_col=0)
    input_x = input_x.iloc[:, 0:sum(config['model']['x_in_dim'])]
    input_x_var = pd.read_csv(f'{data_dir}/input.x_var.csv', index_col=0)
    for index in range(config['model']['x_in_dim'][0]):
        input_x.iloc[:, index] = input_x.iloc[:, index] / input_x_var.iloc[:, index]
    input_x_label = pd.read_csv(f'{data_dir}/input.label.csv', index_col=0)
    train_dataset = DataSet_v3(x_in=input_x,
                               x_out=input_x_out,
                               x_var=input_x_var,
                               label=input_x_label,
                               x_dim=config['model']['x_dim'],
                               warmup_flag=True)
    train_dataloader = torch.utils.data.DataLoader(train_dataset, batch_size=config['train']['batch_size'])
    try:
        torch.manual_seed(config['train']['seed'])
    except KeyError:
        torch.manual_seed(1)
    try:
        gamma_alpha_kl_param=config['model']['gamma_alpha_kl_param']
    except KeyError:
        gamma_alpha_kl_param=None
    try:
        gamma_beta_kl_param=config['model']['gamma_beta_kl_param']
    except KeyError:
        gamma_beta_kl_param=None
    model = model_name(x_in_dim=config['model']['x_in_dim'],
                       x_out_dim=config['model']['x_out_dim'],
                       x_out_likelihood=config['model']['x_out_likelihood'],
                       x_out_mixture=config['model']['x_out_mixture'],
                       h_in_shared_dim=config['model']['h_in_shared_dim'],
                       h_in_z_dim=config['model']['h_in_z_dim'],
                       h_in_y_dim=config['model']['h_in_y_dim'],
                       z_dim=config['model']['z_dim'],
                       y_dim=config['model']['y_dim'],
                       y_kl_param=config['model']['y_kl_param'],
                       y_kl_param_learnable=config['model']['y_kl_learnable'],
                       gamma_kl_param_learnable=config['model']['gamma_kl_learnable'],
                       gamma_fixed=config['model']['gamma_fixed'],
                       gamma_alpha_kl_param=gamma_alpha_kl_param,
                       gamma_beta_kl_param=gamma_beta_kl_param,
                       h_out_shared_dim=config['model']['h_out_shared_dim'],
                       h_out_dims=config['model']['h_out_dims'],
                       batch_norm=config['model']['batch_norm'],
                       layer_norm=config['model']['layer_norm'])

    writer = SummaryWriter(f'{model_dir}/Log/')
    debug_dir = f'{model_dir}/Debug/'
    os.makedirs(debug_dir, exist_ok=True)
    # First run warm-up steps
    optimizer = torch.optim.Adam(list(model.parameters()), amsgrad=True, lr=config['train']['warmup']['learning_rate'])
    likelihood_weight = np.array(config['train']['warmup']['likelihood_weight'])
    likelihood_weight = likelihood_weight / sum(model.x_in_dim)
    likelihood_weight = torch.tensor(likelihood_weight,
                                     dtype=torch.float)
    model.to('cuda')
    counter = 0
    epoch_losses = []
    for i in range(config['train']['warmup']['steps']):
        start = time.time()
        epoch_loss, counter = run_one_epoch(counter, True, train_dataloader, model, likelihood_weight,
                                            optimizer, writer, debug_dir,
                                            beta=repeat(config['train']['warmup']['beta']),
                                            epsilon=repeat(config['train']['warmup']['epsilon'])
                                            )
        end = time.time()
        print(f"Warmup Epoch {i} loss: {epoch_loss}, time elapsed: {end-start}")
        epoch_losses.append(epoch_loss)
    torch.save(model.state_dict(), f'{model_dir}/model.after.warmup.pt')
    # next run decay steps
    train_dataset = DataSet_v3(x_in=input_x,
                               x_out=input_x_out,
                               x_var=input_x_var,
                               label=input_x_label,
                               x_dim=config['model']['x_dim'],
                               warmup_flag=False)
    train_dataloader = torch.utils.data.DataLoader(train_dataset, batch_size=config['train']['batch_size'])

    optimizer = torch.optim.Adam(list(model.parameters()), amsgrad=True,
                                 lr=config['train']['decay']['initial_learning_rate'])
    lr_scheduler = PolynomialLRDecay(optimizer=optimizer,
                                     max_decay_steps=config['train']['decay']['max_decay_steps'],
                                     end_learning_rate=config['train']['decay']['end_learning_rate'])
    likelihood_weight = np.array(config['train']['decay']['likelihood_weight'])
    likelihood_weight = likelihood_weight / sum(model.x_in_dim)
    likelihood_weight = torch.tensor(likelihood_weight,
                                     dtype=torch.float)
    for i in range(config['train']['decay']['steps']):
        start = time.time()
        epoch_loss, counter = run_one_epoch(counter, True, train_dataloader, model,
                                            likelihood_weight, optimizer, writer, debug_dir,
                                            beta=repeat(config['train']['decay']['beta']),
                                            epsilon=repeat(config['train']['decay']['epsilon'])
                                            )
        end = time.time()
        print(f"Decay Epoch {i} loss: {epoch_loss}, time elapsed: {end-start}")
        epoch_losses.append(epoch_loss)
        lr_scheduler.step()
    torch.save(model.state_dict(), f'{model_dir}/model.pt')
    plt.figure()
    plt.plot(epoch_losses)
    plt.savefig(f'{model_dir}/epoch.loss.pdf')
    return model


def test_model(config, counter=None):
    model_name = DeepGenerativeMixModel
    data_dir = config['input']['data_dir']
    model_dir = config['output']['model_dir']
    input_x = pd.read_csv(f'{data_dir}/input.x.csv', index_col=0)
    input_x = input_x.iloc[:, 0:sum(config['model']['x_in_dim'])]
    input_x_out = pd.read_csv(f'{data_dir}/input.x.csv', index_col=0)
    input_x_var = pd.read_csv(f'{data_dir}/input.x_var.csv', index_col=0)
    for index in range(config['model']['x_in_dim'][0]):
        input_x.iloc[:, index] = input_x.iloc[:, index] / input_x_var.iloc[:, index]
    input_x_label = pd.read_csv(f'{data_dir}/input.label.csv', index_col=0)

    test_dataset = DataSet_v3(x_in=input_x,
                              x_out=input_x_out,
                              x_var=input_x_var,
                              label=input_x_label,
                              x_dim=config['model']['x_dim'],
                              train_flag=False)
    test_dataloader = torch.utils.data.DataLoader(test_dataset, batch_size=100)

    model = model_name(x_in_dim=config['model']['x_in_dim'],
                       x_out_dim=config['model']['x_out_dim'],
                       x_out_likelihood=config['model']['x_out_likelihood'],
                       x_out_mixture=config['model']['x_out_mixture'],
                       h_in_shared_dim=config['model']['h_in_shared_dim'],
                       h_in_z_dim=config['model']['h_in_z_dim'],
                       h_in_y_dim=config['model']['h_in_y_dim'],
                       z_dim=config['model']['z_dim'],
                       y_dim=config['model']['y_dim'],
                       y_kl_param=config['model']['y_kl_param'],
                       y_kl_param_learnable=config['model']['y_kl_learnable'],
                       gamma_kl_param_learnable=config['model']['gamma_kl_learnable'],
                       gamma_fixed=config['model']['gamma_fixed'],
                       h_out_shared_dim=config['model']['h_out_shared_dim'],
                       h_out_dims=config['model']['h_out_dims'],
                       batch_norm=config['model']['batch_norm'],
                       layer_norm=config['model']['layer_norm'])

    if counter is None:
        model.load_state_dict(torch.load(f'{model_dir}/model.pt'))
    else:
        model.load_state_dict(torch.load(f'{model_dir}/Debug/model.{counter}.pt'))
    model.to('cuda')
    model.eval()

    likelihood_weight = np.array(config['train']['decay']['likelihood_weight'])
    likelihood_weight = likelihood_weight / sum(model.x_in_dim)
    likelihood_weight = torch.tensor(likelihood_weight,
                                     dtype=torch.float)

    logBFs, y_logits, z_mus, z_sigmas, reconstruction_means, reconstruction_vars = output_log_BF(test_dataloader, model, likelihood_weight)
    logBFs = pd.DataFrame(logBFs, index=input_x.index)
    reconstruction_means = pd.DataFrame(reconstruction_means, index=input_x.index)
    reconstruction_vars = pd.DataFrame(reconstruction_vars, index=input_x.index)

    y_logits = pd.DataFrame(y_logits, index=input_x.index)
    z_mus = pd.DataFrame(z_mus, index=input_x.index)
    z_sigmas = pd.DataFrame(z_sigmas, index=input_x.index)

    if counter is None:
        logBFs.to_csv(f'{model_dir}/log_BFs.csv')
        reconstruction_means.to_csv(f'{model_dir}/reconstruction_means.csv')
        reconstruction_vars.to_csv(f'{model_dir}/reconstruction_vars.csv')
        y_logits.to_csv(f'{model_dir}/y_logits.csv')
        z_mus.to_csv(f'{model_dir}/z_mus.csv')
        z_sigmas.to_csv(f'{model_dir}/z_sigmas.csv')
    else:
        logBFs.to_csv(f'{model_dir}/log_BFs.{counter}.csv')
        reconstruction_means.to_csv(f'{model_dir}/reconstruction_means.{counter}.csv')
        reconstruction_vars.to_csv(f'{model_dir}/reconstruction_vars.{counter}.csv')
        y_logits.to_csv(f'{model_dir}/y_logits.{counter}.csv')
        z_mus.to_csv(f'{model_dir}/z_mus.{counter}.csv')
        z_sigmas.to_csv(f'{model_dir}/z_sigmas.{counter}.csv')
        return logBFs, y_logits, z_mus, z_sigmas, reconstruction_means, reconstruction_vars


def continue_train_model(data_dir, model_dir, counter=500, n_epoch=500):
    input_x = pd.read_csv(f'{data_dir}/input.x.csv', index_col=0)
    input_x_var = pd.read_csv(f'{data_dir}/input.x_var.csv', index_col=0)
    input_x_label = pd.read_csv(f'{data_dir}/input.label.csv', index_col=0)
    train_dataset = DataSet(x=input_x,
                            x_var=input_x_var,
                            label=input_x_label,
                            x_dim=[80, 2])
    train_dataloader = torch.utils.data.DataLoader(train_dataset, batch_size=2000)
    torch.manual_seed(0)

    model = DeepGenerativeMixModel(x_in_dim=[80, 2],
                                   x_out_likelihood=['Binomial', 'Poisson']
                                   )
    model.load_state_dict(torch.load(f'{model_dir}/model.pt'))
    model.to('cuda')

    writer = SummaryWriter(f'{model_dir}/Log/')
    debug_dir = f'{model_dir}/Debug/'
    os.makedirs(debug_dir, exist_ok=True)
    optimizer = torch.optim.Adam(list(model.parameters()), amsgrad=True, lr=0.0001)
    model.to('cuda')
    epoch_losses = []
    for i in range(n_epoch):
        start = time.time()
        epoch_loss, counter = run_one_epoch(counter, True, train_dataloader, model, optimizer, writer, debug_dir)
        end = time.time()
        print(f"Epoch {i} loss: {epoch_loss}, time elapsed: {end-start}")
        epoch_losses.append(epoch_loss)
    torch.save(model.state_dict(), f'{model_dir}/model.continue.pt')
    plt.figure()
    plt.plot(epoch_losses)
    plt.savefig(f'{model_dir}/epoch.loss.continue.pdf')
    return model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    parser.add_argument('--mode', type=int, default=0)
    parser.add_argument('--random', type=str, default='0')
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    if args.mode == 0:
        train_model(config)
    elif args.mode == 1:
        test_model(config)


if __name__ == '__main__':
    main()
