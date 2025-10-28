# Dengue vaccine comparison 

This repository contains the code needed to reproduce the main results presented in the manuscript: Modelling the efficacy of Qdenga, Dengvaxia, and the Butantan-DV dengue vaccines: a comparative analysis and open questions<img width="468" height="127" alt="image" src="https://github.com/user-attachments/assets/94dde4da-bc13-4f40-a845-184127adabd8" />
.

Instructions to download and install cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan). Instructions to download and install Rtools can be found [here](https://cran.r-project.org/bin/windows/Rtools/).

## Repository structure


### Scripts
* *1.fit_all_trials.R*: # Script to calibrate the stan model to the publicly available data of all three vaccines. The script is set up to run the final model (M21) first and all other model variants presented in the manuscript can be run  below. All model variants are run using the same stan model (final_model.stan) with different parameters turned on and off, using the flags. 
* *2_plot_VE_figure.R*: Main script to plot the attack rates and vaccine efficacy estimated in *1_VE_model_fitting.R*, to reconstruct Figure 1 in the manuscript. 

### Data

* *hosp_data.csv*: hospital VCD case data extracted from the published phase III clinical trial [1-5]. 
* *vcd_data.csv*: symptomatic VCD case extracted from the published phase III clinical trial [1-5]. 
* *n0_new.csv*: initial neutralising antibody titres induced by Qdenga vaccination, fitted from the published phase III clinical trial [1-5]. 
* *seropositive_by_age_baseline.csv*: number of baseline seropositive individuals by age group, extracted from the published phase III clinical trial [1-5].

### Models

* *final_model.csv*: Rstan model to fit the Bayesian survival model to case data, assuming binomial and multinomial likelihoods. Flags in the data block turn on and off different parameters used to test the 31 model variants presented in the manuscript. 


### R
* *factor_data.R*: functions to format the data.
* *format_stan_data.R*: function to format all data input to the Stan model used for model calibration. Outputs data in a list format.
* *run_model_fitting.R*: function to calibrate the stan model to the publically available data.
* *plot_model_outputs.R*: functions to calculate the observed attack rates and format model outputs for plotting. 
