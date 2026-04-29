# ------------------------------------------------------------------------------
# extract_model_results
#
# Convert Stan posterior draws into summarised model estimates by computing
# posterior means and credible intervals for each parameter. 
# ------------------------------------------------------------------------------

extract_model_results = function(fit_ext){
  out =  fit_ext %>%  
    as.data.frame() %>%  
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value * 100) %>%  
    summarise(
      lower = quantile(value, 0.025, na.rm = T),
      mean = mean(value, na.rm = T),
      upper = quantile(value, 0.975, na.rm = T)
    )
  
  out2 = out %>%  
    separate(name, into = c("group", "name"), sep = "\\[") %>% 
    separate(name, into = c("name", NA), sep = "\\]") %>% 
    mutate(type = "model") 
  
  return(out2)
}


# ------------------------------------------------------------------------------
# process_VE_files
#
# Read and foramt files across sensitivity analyses
# ------------------------------------------------------------------------------

process_VE_files = function(folder_names, type = c("Q", "De", "Bu"), n) {
  type = match.arg(type)
  
  # build file paths
  files = file.path("output", folder_names, paste0("VE_", type, ".RDS"))
  
  # check existence
  exists = file.exists(files)
  files = files[exists]
  
  # read files
  VE_list = lapply(files, readRDS)
  names(VE_list) = folder_names[exists]
  
  # optionally select only VE columns (skip for B if not needed)
  if (type %in% c("Q")) {
    VE_list = lapply(VE_list, function(df) {
      dplyr::select(df, dplyr::contains(paste0("VE_", type)))
    })
  }
  
  if (type %in% c("De")) {
    VE_list = lapply(VE_list, function(df) {
      dplyr::select(df, dplyr::contains("VE_D"))
    })
  }
  
  # sample from posterior
  VE_list = lapply(VE_list, function(df) {
    as.data.frame(sapply(df, sample, n))
  })
  
  return(VE_list)
}

# ------------------------------------------------------------------------------
# process_sens
#
# Convert Stan posterior draws across sensitivity analyses into tidy format 
# ------------------------------------------------------------------------------

process_sens = function(list, split_into) {
  list %>% 
    dplyr::bind_rows(.id = "sa") %>% 
    dplyr::mutate(ni = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -c(ni, sa)) %>%
    dplyr::group_by(name) %>% 
    dplyr::mutate(value = value * 100) %>%
    tidyr::separate(name, into = c("group", "name"), sep = "\\[") %>%
    tidyr::separate(name, into = c("name", NA), sep = "\\]") %>%
    tidyr::separate(name, into = split_into)
}



# ------------------------------------------------------------------------------
# calculate_diff
#
# Calculate the absolute difference in age-weighted mean VE between
# each SA model and the baseline model 
# ------------------------------------------------------------------------------


calculate_diff = function(df,
                           age_weights = NULL,
                           main_model,
                           outcome = 1,
                           serostatus_labels) {

  # filter to baseline scenario
  baseline = df %>%
    filter(outcome == outcome, sa == main_model)
  
  # take weighted average across age groups 
  baseline_mean = baseline %>%
    left_join(age_weights, by = "age") %>%
    group_by(serostatus, serotype, month, ni) %>%
    summarise(value = weighted.mean(value, w = n), .groups = "drop") %>%
    group_by(serostatus, serotype) %>%
    summarise(baseline_mean = mean(value), .groups = "drop")
  
  
  # sens means
  scenario = df %>%
    filter(outcome == outcome) %>%
    left_join(age_weights, by = "age") %>%
    group_by(sa, serostatus, serotype, month, ni) %>%
    summarise(value = weighted.mean(value, w = n), .groups = "drop") %>%
    group_by(sa, serostatus, serotype) %>%
    summarise(scenario_mean = mean(value), .groups = "drop")
  
  # join and compute difference
  scenario %>%
    left_join(baseline_mean, by = c("serostatus", "serotype")) %>%
    mutate(
      difference  = (scenario_mean - baseline_mean),
      serotype   = factor(paste0("DENV", serotype)),
      serostatus = factor(
        serostatus, levels = 1:3,
        labels = unname(serostatus_labels)
      )
    )  %>%  # include baseline for reference
  mutate(difference = ifelse(sa == main_model, baseline_mean, difference),
         main = ifelse(sa == main_model, T, F))  
    
}



# ------------------------------------------------------------------------------
# plot_heatmap
#
# Plot heatmap of vaccine efficacy differences between sensitivity analyses and
# the baseline model, stratified by serostatus, serotype, and clinical outcome.
# ------------------------------------------------------------------------------

plot_heatmap = function(data, scale_limit) {
  

  ggplot(data, aes(x = serotype, y = sa, fill = difference)) +
    geom_tile(
      data = filter(data, main),
      fill = "darkgrey",
      colour = "white",
      linewidth = 0.4
    ) +
    geom_tile(
      data = filter(data, !main),
      colour = "white",
      linewidth = 0.4
    ) +
    geom_text(aes(
      label    = sprintf("%+.1f", difference),
      fontface = ifelse(main, "bold", "plain")
    ), size = 2.5) +
    scale_fill_gradient2(
      low      = "#185FA5",
      mid      = "grey97",
      high     = "#C44E52",
      midpoint = 0,
      limits   = c(-scale_limit, scale_limit),
      breaks   = seq(-scale_limit, scale_limit, by = 10),
      name     = "Absolute\ndifference from\nbaseline (%)"
    ) +
    facet_wrap( ~ serostatus, nrow = 1, ncol = 3) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    my_theme +
    scale_y_discrete(limits = rev) 
  
}

