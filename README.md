# Dengue vaccine comparison 

This repository contains the code needed to reproduce the main results presented in the manuscript: Modelling the efficacy of Qdenga, Dengvaxia, and the Butantan-DV dengue vaccines: a comparative analysis and open questions.

Instructions to download and install cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan). Instructions to download and install Rtools can be found [here](https://cran.r-project.org/bin/windows/Rtools/).

## Repository structure


### Scripts
* *1.fit_all_trials.R*: Script to calibrate the stan model to the publicly available data of all three vaccines. The script is set up to run the final model (M21) first and all other model variants presented in the manuscript can be run  below. All model variants are run using the same stan model (final_model.stan) with different parameters turned on and off, using the flags. 
* *2.plot_all_attack_rates.R*: Script to plot the final model fits of all three vaccines. Script outputs Fig.1.
* *3.plot_efficacy.R* Script to plot 24 month VE for all vaccines aggregated over time. Script outputs Fig.2
* *4.plot_efficacy_over_time* Script to plot VE of Dengvaxia and Qdenga up to 54 months. Script outputs Figs 3 and 4. 
  
### Data

### Qdenga 
* *processed_Q/hosp_data.csv*: hospital VCD cases extracted from the published phase III clinical trial. 
* *processed_Q/vcd_data.csv*: symptomatic VCD cases extracted from the published phase III clinical trial.
* *processed_Q/n0_new.csv*: initial neutralising antibody titres induced by Qdenga vaccination, fitted from the published phase III clinical trial.
* *processed_Q/seropositive_by_age_baseline.csv*: number of baseline seropositive individuals by age group, extracted from the published phase III clinical trial.

#### Dengvaxia
* *processed_De/cases_stan_format.RDS*: symptomatic, hospitalised, severe VCD cases extracted from the published phase III clinical trial.
* *processed_De/mu.csv*: initial neutralising antibody titres induced by Dengvaxia vaccination, fitted from the published phase III clinical trial.
* *processed_De/baseline_SP.csv*: number of baseline seropositive individuals by age group, extracted from the published phase III clinical trial.

#### Butantan-DV

* *processed_Bu/cases_stan_format.RDS*: symptomatic VCD cases extracted from the published phase III clinical trial.
* *processed_Bu/mu.RDS*: initial neutralising antibody titres induced by Butantan-DV vaccination, from the published phase III clinical trial.
* *processed_Bu/baseline_SP.RDS*: number of baseline seropositive individuals by age group, extracted from the published phase III clinical trial.


### Models

* *final_model.stan*: Rstan model to fit the Bayesian survival model to case data, assuming binomial and multinomial likelihoods. Flags in the data block turn on and off different parameters used to test the 21 model variants presented in the manuscript.
* *final_model_severe.stan*: *final_model.stan* but additionally fits Dengvaxia severe data. 

### R
* *factor_data_Q.R*: functions to format the Qdenga data.
* *format_stan_data.R*: vaccine-specific functions to format all data input to the Stan model used for model calibration. Outputs data in a list format.
* *plot_model_outputs.R*: vaccine-specific functions to calculate the observed attack rates and format model outputs for plotting. 
* *plotting_functions.R*: functions to plot all vaccine efficacies together. 
* *run_model_fitting_joint.R*: function to calibrate the stan model to the publically available data.
