# Dengue vaccine comparison

This repository contains the code needed to reproduce the main results presented in the manuscript:

**_Modelling the efficacy of Qdenga, Dengvaxia, and the Butantan-DV dengue vaccines: a comparative analysis and open questions._**

The workflow combines publicly available clinical trial datasets with a Bayesian survival model implemented in Stan. Trial datasets are standardized, reshaped into fixed-dimension arrays, and merged into a unified Stan input structure before model fitting and posterior analysis.

Instructions to download and install **CmdStan** can be found [here](https://mc-stan.org/users/interfaces/cmdstan).  
Instructions to download and install **Rtools** (Windows only) can be found [here](https://cran.r-project.org/bin/windows/Rtools/).


---

## Overview of workflow

The pipeline follows a reproducible sequence:

1. **Import & standardise trial data**  
   Raw datasets are converted to consistent categorical formats to ensure reproducible aggregation and plotting.

2. **Trial-specific data formatting**  
   Each vaccine dataset is summarised by age, serostatus, serotype, trial arm, and time. Counts are reshaped into multidimensional arrays required by Stan. Vaccine-specific formatted datasets are merged into a single Stan data list. 

3. **Model calibration**  
   The Stan model is compiled and sampled using `cmdstanr`. Posterior draws, derived quantities, diagnostics, and evaluation metrics are saved.

4. **Posterior analysis & visualisation**  
   Model outputs are converted into interpretable vaccine efficacy and attack rate summaries for plotting.

Scripts should be run in order (1–5) to reproduce manuscript figures.


## Repository structure

### Scripts

Scripts execute the full pipeline from model calibration to figure generation:

- **1.fit_all_trials.R**  
  Runs the multi-trial Bayesian calibration using all vaccine datasets. The script is configured to run the final model (M21) first. Alternative model variants are implemented via configuration flags passed to the Stan model.

- **2.plot_all_attack_rates_scatter.R**  
  For each vaccine, joins posterior attack rate estimates to observed trial data and produces observed vs. model-estimated scatter plots. Saves a combined faceted ggplot object (by vaccine) to disk for use in the next script.  
  → Outputs `output/figures/all_attack_rates.rds`

- **3.plot_all_attack_rates.R**  
  Loads the scatter plot produced by script 2 and combines it with per-vaccine time-course and serotype attack rate plots into a single assembled figure.  
  → Outputs manuscript **Fig. 1** (`output/figures/Fig1.pdf`)

- **4.plot_efficacy_over_time.R**  
  Plots model-estimated vaccine efficacy trajectories up to 54 months post final dose, stratified by serostatus, serotype, and age group, for all three vaccines. Produces separate figures for symptomatic and hospitalised outcomes.  
  → Outputs manuscript **Fig. 2** (`output/figures/Fig2.pdf`) and **Fig. 3** (`output/figures/Fig3.pdf`)

- **5.calculate_plot_ve_sensitivity_analysis.R**  
  Loads posterior VE draws from the baseline model (M4) and a set of structural sensitivity analysis variants. Computes age-weighted VE differences between each alternative model and the baseline, and visualises results as diverging heatmaps stratified by serostatus, serotype, and outcome. Individual vaccine heatmaps and a combined figure are saved.  
  → Outputs manuscript **Fig. 4** (`output/figures/Fig4.pdf`) and individual heatmaps per vaccine


### Data

All datasets are preprocessed aggregations of published phase III clinical trial results. 

Where possible case counts are extracted by month, serotype, infecting serotype, age group and clinical outcome (symptomatic, hospitalised, and severe VCD).

#### Qdenga

- `processed_Q/hosp_data.csv` — hospitalised VCD cases  
- `processed_Q/vcd_data.csv` — symptomatic VCD cases  
- `processed_Q/n0_new.csv` — initial fitted neutralising titres  
- `processed_Q/seropositive_by_age_baseline.csv` — age-stratified baseline serostatus counts  

#### Dengvaxia

- `processed_De/cases_stan_format.RDS` — symptomatic/hospital/severe VCD cases
- `processed_De/mu.csv` — initial fitted neutralising titres  
- `processed_De/baseline_SP.csv` — baseline serostatus counts  

#### Butantan-DV

- `processed_Bu/cases_stan_format.RDS` — symptomatic VCD cases  
- `processed_Bu/mu.RDS` — neutralizing titres  
- `processed_Bu/baseline_SP.RDS` — baseline serostatus counts  

---

### Models

Stan models implement the Bayesian survival framework used for calibration.

- **final_model.stan**  
  Core multi-trial survival model with binomial and multinomial likelihoods. Configuration flags enable all model variants tested in the manuscript through script **1.fit_all_trials.R** 

- **final_model_severe.stan**  
  Extension of the main model including Dengvaxia severe outcome likelihoods.

---

### R 


Functions implement data preprocessing, formatting, model fitting, and plotting:

- **format_stan_data.R**  
  Aggregates case counts, reshapes multidimensional arrays, and combines priors and flags into Stan-ready inputs.

- **run_model_fitting_joint.R**  
  High-level wrapper that preprocesses data, constructs Stan inputs, runs MCMC sampling via `cmdstanr`, and saves posterior outputs and diagnostics.

- **format_model_outputs_for_plotting.R**  
  Formats model outputs (attack rates or vaccine efficacy) ready for plotting. Calculates absolute difference in vaccine efficacy estimates for sensitivity analysis. 

- **format_data_for_plotting.R**  
  Formats raw case data into observed attack rates ready for plotting.
  
- **calculate_bayes.R**
  Computes Savage-Dickey Bayes factors for the waning parameter (L) across all three vaccines. Estimates prior and posterior densities at zero using a truncated normal prior and logspline density estimation, and writes the results to BF.csv.
 

---

## Output

Model fitting generates:

- posterior chains and summaries  
- vaccine efficacy and attack rate posterior estimates  
- WAIC/log-likelihood diagnostics  

Figure generation produces PDFs in `output/figures/`:

| File | Script | Description |
|------|--------|-------------|
| `all_attack_rates.rds` | 2 | Intermediate: combined scatter plot object |
| `Fig1.pdf` | 3 | Observed vs. estimated attack rates |
| `Fig2.pdf` | 4 | VE against symptomatic disease over time |
| `Fig3.pdf` | 4 | VE against hospitalisation over time |
| `Fig4.pdf` | 5 | Sensitivity analysis heatmaps |
| `Q_plot_heatmap_diff.pdf` | 5 | Qdenga sensitivity heatmap (standalone) |
| `D_plot_heatmap_diff.pdf` | 5 | Dengvaxia sensitivity heatmap (standalone) |
| `B_plot_heatmap_diff.pdf` | 5 | Butantan-DV sensitivity heatmap (standalone) |

