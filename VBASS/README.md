## VBASS
  
  VBASS is a new model to integrate single cell genomics data in identifying de novo risk genes. 

## Environment:
  
  `conda env create -f MLenv.yml`
  
  `conda activate MLenv`

## Usage:

  Prepare the desired input, see `prepare.input.folders.R/` for examples.
  
  Prepare the running configs into a `.json` file, see `config.json/real/config.real.json` for examples.
  
  Run `python train_mix_model.py --config config.json` to train the model. It takes ~20 mins on a Nvidia A40 GPU.
  
  Run `python train_mix_model.py --config config.json --mode 1` to output the results.
  
## Result:
  
  A folder with several files:
  
  `y_logits.csv`, the second column is the estimated π in log scale.
  
  `log_BFs.csv`, the Bayes factor, or log ratio of posterior probability of alternative hypothesis over null hypothesis in log scale.
  
  `reconstruction_vars.csv`, estimated effect sizes variance.
  
  `reconstruction_means.csv`, estimated effect sizes mean.
  
  `z_sigmas.csv, z_mus.csv`, latent representation of single cells.
  
  `epoch.loss.pdf`, training loss curve.
  
  
## Downstream analysis:

  see `fig.*.R` for examples of downstream analysis to reproduce the figures in our manuscript.
  
