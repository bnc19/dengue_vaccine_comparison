# function to add extra populations to VCD data and calculate AR --------------- 

calc_CYD_attack_rates = function(VCD, titre = F) {
  out = VCD %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  
  return(out)
}

# function to extract model results --------------------------------------------
extract_CYD_model_results = function(fit_ext){
  
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


# function to plot attack rates ------------------------------------------------

plot_CYD_attack_rate = function(VCD,
                            file_path, 
                            hospital = F,
                            AR) {
  
# add aggregated populations to data and calculate attack rates
AR_data = lapply(VCD, calc_CYD_attack_rates)
AR_model = extract_CYD_model_results(AR)  

age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = scales::brewer_pal(palette = "RdPu")(6)[2:5]
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]

# plot symp attack rate by serotype and trial arm ------------------------------
  
  Sy_AR_plot_VK = AR_model %>%
    filter(group == "V_AR_VK") %>%
    separate(name, into = c("arm", "serotype")) %>% 
    mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
           serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4"))) %>% 
    bind_rows(AR_data$Sy_VK) %>%
    ggplot(aes(x = arm, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "") +
    scale_color_manual(values = serotype_fill) +
   guides(shape = "none", linetype = "none")
  
# plot symp attack rate by age, serostatus and trial arm -----------------------
  Sy_AR_plot_BVJ = AR_model %>%
    filter(group == "V_AR_BVJ") %>%
    separate(name, into = c("serostatus", "arm", "age")) %>% 
    mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
           arm = factor(arm, labels = c("placebo", "vaccine")),
           age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
    bind_rows(AR_data$Sy_BVJ) %>%
    unite(c(arm, serostatus), col = "x", sep = "\n") %>%
    ggplot(aes(x = x, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = age,
        group = interaction(type, age)),
      position = position_dodge(width = 0.5),
      size = 3) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, age),
        linetype = type,
        color = age),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1) +
    labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  scale_color_manual(values = age_fill)
  
# plot hosp attack rate by age, serostatus, serotype and trial arm -------------

  H_AR_plot_BVKJ = AR_model %>%
  filter(group == "H_AR_BVKJ") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
           serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4")),
           serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
           age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Ho_BVKJ) %>%
  unite(c(arm, serostatus), col = "x", sep = "\n") %>%
  ggplot(aes(x = x, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = serotype_fill) +
    theme(legend.position = "none") +
    facet_wrap(~ age) 
  
# plot hosp attack rate by age, serostatus, trial arm and time -----------------

  H_AR_plot_BVJD =  AR_model %>%
      filter(group == "H_AR_BVJD") %>%
      separate(name, into = c("serostatus", "arm", "age", "time")) %>% 
      mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
             arm = factor(arm, labels = c("placebo", "vaccine")),
             age = factor(age, labels = c("2-8yrs", "9-16yrs")),
             time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60"))) %>%
      bind_rows(AR_data$Ho_BVJD ) %>%
      ggplot(aes(x = time, y = mean)) +
      geom_point(
        aes(
          shape = type,
          color = age,
          group = interaction(type, age)
        ),
        position = position_dodge(width = 0.7),
        size = 3
      ) +
      geom_errorbar(
        aes(
          ymin = lower ,
          ymax = upper ,
          group = interaction(type, age),
          linetype = type,
          color = age
        ),
        position = position_dodge(width =  0.7),
        width =  0.4,
        linewidth = 1
      ) +
      facet_grid(serostatus ~ arm) +
      labs(x = "Month", y = "Hospitalisation \nattack rate (%)") +
      scale_color_manual(values = age_fill) +
      theme(legend.position = "none")  
    
  out = cowplot::plot_grid(
    cowplot::plot_grid( 
      Sy_AR_plot_BVJ,
      Sy_AR_plot_VK, ncol =2,
      labels = c("a", "b")), 
    H_AR_plot_BVJD,
    H_AR_plot_BVKJ,
    ncol = 1,
    rel_heights = c(1,1.3,1), 
    labels = c("", "c", "d")
  )
  
    ggsave(
      plot = out,
      filename = paste0(file_path, "/AR.png"),
      height = 60,
      width = 70,
      units = "cm",
      dpi = 300,
      scale = 0.75
    )
    
    return(out)
}

# function to plot severe AR  --------------------------------------------------
plot_severe_attack_rate = function(VCD,
                                   file_path, 
                                   AR) {
  
  # add aggregated populations to data and calculate attack rates
  sev = filter(VCD, outcome == "Se")
  
  AR_data = calc_CYD_attack_rates(sev)
  AR_model = extract_CYD_model_results(AR)  
  
  age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
  serotype_fill = scales::brewer_pal(palette = "RdPu")(6)[2:5]
  trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]
  
  # plot severe attack rate by serotype, trial arm, serotype and age
  
  Se_AR_BVKJ = AR_data %>% 
    filter(time == "all") %>%  
    bind_rows(AR_model %>%
                filter(group == "S_AR_BVKJ") %>%
                separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
                mutate(arm = factor(arm, labels=c("placebo", "vaccine")),
                       serotype = factor(serotype, labels =c("DENV1","DENV2","DENV3", "DENV4")),
                       serostatus = factor(serostatus, labels = c("seronegative", "seropositive")), 
                       age = factor(age, labels = c("2-8yrs", "9-16yrs")))) %>%  
    filter(!(arm == "vaccine" & serostatus == "seropositive")) %>%  
    ggplot(aes(x = arm, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Severe attack rate (%)") +
    scale_color_manual(values = serotype_fill) +
    guides(shape = "none",
           linetype = "none") 
    facet_grid(age ~ serostatus)
  
  # plot severe attack rate by age, arm, time in SN   ----------------------------
  
  Se_AR_VJD = AR_data %>% 
    filter(outcome == "Se", time != "all", serostatus == "seronegative") %>%  
    bind_rows(AR_model %>%
                filter(group == "S_AR_snVJD") %>%
                separate(name, into = c("arm", "age", "time")) %>% 
                mutate(arm = factor(arm, labels=c("placebo", "vaccine")),
                       time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60")), 
                       age = factor(age, labels = c("2-8yrs", "9-16yrs")))) %>%  
    ggplot(aes(x = time, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = arm,
        group = interaction(type, arm)
      ),
      position = position_dodge(width = 1),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, arm),
        linetype = type,
        color = arm
      ),
      position = position_dodge(width =  1),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Severe attack rate (%)") +
    guides(shape = "none",
           linetype = "none") +
    facet_wrap(~age) +
    scale_color_manual(values = trial_fill)
  
  
# plot severe attack rate by age in vac SP -------------------------------------
  
Se_AR_SP_V = AR_data %>% 
    filter(outcome == "Se", arm == "vaccine", serostatus == "seropositive") %>%  
    bind_rows(AR_model %>%
                filter(group == "S_AR_spvJ") %>%
                separate(name, into = c("age")) %>% 
                mutate(age = factor(age, labels = c("2-8yrs", "9-16yrs")))) %>%  
    ggplot(aes(x = age, y = mean)) +
    geom_point(
      aes(shape = type),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        linetype = type
      ),
      position = position_dodge(width = 0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Severe attack rate (%)") +
    guides(shape = "none",
           linetype = "none")
  
  ggsave(
    plot = cowplot::plot_grid(
      Se_AR_VJD,Se_AR_SP_V, 
      Se_AR_BVKJ,labels = c("a", "b", "c"), ncol = 1), 
    filename = paste0(file_path, "/severe_AR.png"),
    height = 40,
    width = 50,
    units = "cm",
    dpi = 300,
    scale = 0.75
  )
  
}

# function to plot vaccine efficacy --------------------------------------------

plot_CYD_VE  = function(file_path,
                    hospital = F,
                    VE = NULL) {
  
  VE_model = extract_CYD_model_results(VE)  
  
  VE_model$lower = ifelse(VE_model$lower < -400, -400, VE_model$lower)
  if(hospital ==F){ 
    fil = "1"
    y = "Efficacy against symptomatic disease (%)"
  } else {
    fil= "2" 
    y = "Efficacy against hospitalisation (%)"
  }
  
  age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
  serotype_fill = scales::brewer_pal(palette = "RdPu")(6)[2:5]
  trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]
  
  # plot VE by serostatus for each serotype --------------------------------------
  
  VE_BKT =  VE_model %>%
    filter(group == "VE_BKRT") %>%
    separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>% 
    filter(outcome == eval(rlang::parse_expr(fil)), serostatus != 3) %>% 
    mutate(serotype = factor(serotype, 
                             labels = c("DENV1","DENV2", "DENV3","DENV4")),
           month = as.numeric(month),
           serostatus = factor(serostatus, 
                               labels = c("seronegative", "monotypic")))  
  
  VE_BKT_plot = VE_BKT %>% 
    ggplot(aes(x = month , y = mean)) +
    geom_line(aes(color = serostatus)) +
    geom_ribbon(
      aes(
        ymin = lower ,
        ymax = upper ,
        fill = serostatus
      ), alpha = 0.5
    ) +
    labs(x = "Month", y = y) +
    scale_x_continuous(breaks = seq(0, 54, 12)) +
    scale_y_continuous(breaks = seq(round(min(VE_BKT$lower),-2), 100, 50)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    facet_grid(serostatus~serotype) + 
    theme_light() + 
    theme(legend.position = "top") +
    scale_color_manual(values =  trial_fill) +
    scale_fill_manual(values = trial_fill) 
  
  ggsave(
    plot = VE_BKT_plot,
    filename = paste0(file_path, "/VE_BKT_plot_outcome_", fil, ".png"),
    height = 21,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.8
  )
  
  # plot VE  age, serostatus for each serotype  ----------------------------------
  
  VE_BKJT =  VE_model %>%
    filter(group == "VE") %>%
    separate(name, into = c("serostatus", "serotype", "age","outcome", "month")) %>%
    filter(outcome == eval(rlang::parse_expr(fil)), serostatus != 3) %>% 
    mutate(
      age = factor(age, labels = c("2-8yrs","9-16yrs")),
      serotype = factor(serotype, 
                        labels = c("DENV1","DENV2", "DENV3","DENV4")), 
      month = as.numeric(month),
      serostatus = factor(serostatus, 
                          labels = c("seronegative", "monotypic")))
  
  
  VE_BKJT_plot = VE_BKJT %>% 
    ggplot(aes(x = month , y = mean)) +
    geom_line(
      aes(color = serostatus)) +
    geom_ribbon(
      aes(
        ymin = lower ,
        ymax = upper ,
        fill = serostatus),
      alpha = 0.5) +
    labs(x = "Month", y = y) +
    scale_x_continuous(breaks = seq(0, 54, 6)) +
    scale_y_continuous(breaks = seq(round(min(VE_BKJT$lower),-2), 100, 50)) +
    geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
    scale_color_manual(values =  trial_fill) +
    scale_fill_manual(values = trial_fill) +
    facet_grid(serotype ~ age) + theme_light() + theme(legend.position = "top")

  
  ggsave(
    plot = VE_BKJT_plot,
    filename = paste0(file_path,"/VE_BKJT_plot_outcome_", fil, ".png"),
    height = 41,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.8
  )
}

# function to plot titres ----------------------------------------------------------------------------------------
plot_CYD_titres = function(file_path, n = NULL){
  
  
  # format model output and add data 
  model_titres = n %>%   
    as.data.frame() %>%  
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = log(value)) %>%  
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) %>%  
    filter(grepl("n\\[", name)) %>% 
    separate(name, into = c(NA, "name"), sep = "\\[") %>% 
    separate(name, into = c("name", NA), sep = "\\]") %>% 
    separate(name, into = c("serostatus", "serotype", "time")) %>%  
    mutate(time= as.numeric(time),
           serostatus = factor(serostatus, labels = c("seronegative", 
                               "seropositive")),
           serotype = factor(serotype, 
                             labels = c("DENV1","DENV2", "DENV3","DENV4")))
  
  
  plot_titres = model_titres %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper,
    ), alpha = 0.2) +
    theme(legend.position ="top") +
    facet_grid(serostatus~ serotype) +
    labs(x = "Month", y = "Log titre")
  
  ggsave(plot_titres, filename = paste0(file_path,"/plot_titres.jpg"),
         width = 40, height = 20, unit ="cm")
}

# function to plot everything --------------------------------------------------
plot_CYD_output = function(file_path,
                           severe = severe) {
  library(tidyverse)
  library(Hmisc)
  library(cowplot)

  # data
  VCD = readRDS("CYD/data/processed/cases_stan_format.RDS")
  AR = readRDS(paste0(file_path, "/AR.RDS"))
  VE = readRDS(paste0(file_path, "/VE.RDS"))
  n = readRDS(paste0(file_path, "/n.RDS"))
  
# plot attack rates
plot_CYD_attack_rate(VCD =  VCD, file_path = file_path, AR = AR)
  
# plot VE
plot_CYD_VE(file_path = file_path, VE = VE)
  
# plot VE against hosp
plot_CYD_VE(file_path = file_path, hospital = T, VE = VE)
  
# plot titres
plot_CYD_titres(file_path = file_path, n = n)
  
# if severe plot severe AR 
if (severe == T)  plot_severe_attack_rate(VCD = VCD, file_path = file_path, AR = AR)

}
