---
editor_options: 
  markdown: 
    wrap: 72
---

# mortality_models

This is a repository for code used to estimate Bayesian models of tree
mortality and survival using periodic Eastwide Forest Inventory Data

## Repository Overview

This repository is ordered into several folders:

1\. `data` - where intial eastwide data is stored

2\. `R code` - R code used to run analyses

3\. `modelcode` - STAN code for Bayesian models

4\. `SPCD_standata_full` - folders for data objects used to run stan
models

5\. `SPCD_stanoutput_full` - folder for outputs from species-level stan
models

6\. `SPCD_stanoutput_joint` - folder for output from hierarchical stan
model

## R code Overview

1.  `R/prepData/`
    -   Folder with Data cleaning Code and data exploration

<!-- -->

2.  `R/speciesModels/`
    -   Folder with code used to estimate Species-level Bayesian models
        for 17 Eastern US species

<!-- -->

3.  `R/speciesBasisModels/`
    -   Folder with code used to estimate Species-level Bayesian models
        with penalized spline on DBH

<!-- -->

4.  `R/hierarchicalModel/`
    -   Folder with code used to estimate the Hierarchical Bayesian
        model for 17 Eastern US species

<!-- -->

5.  `R/modelAssessment/`
    -   Folder with scripts to compare species-level and hierarchical
        model predictive ability, complexity, computational effort, and
        spatial predictions.

## Notes on resources and packages used with this analysis

This repository uses R/Rstudio and `rstan` to sample from several
different Bayesian models. Underlying data is from eastwide Forest
Inventory and Analysis Periodic dataset and is not currently shared on
github. We fit hierarchical model of survival using Cyverse Discovery
Environment (<https://cyverse.org/>).
