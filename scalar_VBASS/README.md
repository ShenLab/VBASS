# scalar_VBASS
scalar_VBASS is new model to integrate functional genomics data in identifying de novo risk genes. In additional to the observation of damage variants, it takes the expression level of genes during embryo development as an input to inform prior of a gene being diesease risk. 

  We hypothesized that disease risk genes should:
  
    1) Harbor enriched damage de novo variants in probands compared to non-affected population. 
    2) Show high expression during corresponding organ development.
    
  The model look like this:
  
  <img src="https://github.com/ShenLab/x-TADA/blob/master/x-TADA.model.png?raw=true" width="512">
  
ϕ_0 and ϕ_1 are hyperparameters to estimate. N is the observed damage variant number in patient cohorts. S is the expression level for each gene, in this CHD practice it is mouse E14.5 heart expression rank percentile.

## Environment:
  
  Install Following packages:
  
  `R>=3.6.0`
  
  `packages: stan, dplyr`
  
  Or use the conda environment we provided:
  
  `conda env create -f rstan.yml`
  
  `conda activate rstan`
  
  Typical install time will be ~10-20 min on a normal desktop computer.

## Usage:

  Source all the R files in `scalar_VBASS/` directory.
  
  Modify the `demo.R` with your input. See input format in demo.R
  
  Run `source(demo.R)`. Expected run time will be ~15 min on a normal desktop computer.
  
## Result:
  
  A list object with three attributes:
  
    1) mcmcDD
  
      mcmcDD raw result
  
    2) pars0
  
      statistical summary of mcmcDD result
  
    3) dataFDR
  
      posterior probability, qvalues, etc for each gene.
      
## Figures in paper
  
  Run `CHD_scalarVBASS_Jinetal.R`. Or you can skip this step and use the results stored in the `result/` folder.
  
  Run `simulation.scalar_VBASS.R <seed>`. Run 100 times with `<seed>` from 1-100. Or you can skip this step and use the results stored in the `RDS.files/` folder.
  
  Run `simulation.extTADA.R <seed>`. Run 100 times with `<seed>` from 1-100. Or you can skip this step and use the results stored in the `RDS.files/` folder.
  
  Run `fig.*.R` files.
  
