# function to add extra populations to VCD data and calculate AR --------------- 
calc_BUT_attack_rates = function(VCD, titre = F) {
  out = VCD %>%  
    mutate(mean =  binconf(Y,N, method ="exact")[,1] * 100,  
           lower = binconf(Y,N, method ="exact")[,2] * 100, 
           upper = binconf(Y,N, method ="exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  
  return(out)
}

# function to extract model results --------------------------------------------
extract_BUT_model_results = function(fit_ext){
  
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

plot_BUT_attack_rate = function(age_cases,
                                serotype_serostatus_cases,
                                file_path,
                                hospital = F,
                                AR) {
  
  
# add aggregated populations to data and calculate attack rates

AR_age_data = calc_BUT_attack_rates(age_cases)
AR_serotype_data = calc_BUT_attack_rates(serotype_serostatus_cases)
AR_model = extract_BUT_model_results(AR)
  
# plot serotype serostatus attack rate 
AR_plot_BVK = AR_model %>%
  filter(group == "AR_BVK") %>%
  separate(name, into = c("Serostatus", "Arm", "Serotype")) %>% 
  mutate(Arm = factor(Arm, levels = 1:2, labels = c("placebo", "vaccine")),
         Serostatus = factor(Serostatus, levels = 1:2, labels = c("seronegative", "seropositive")),
         Serotype = factor(Serotype, levels = 1:2, labels = c(paste0("DENV", 1:2)))) %>% 
    bind_rows(AR_serotype_data) %>%
    ggplot(aes(x = Arm, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = Serotype,
        group = interaction(type, Serotype)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, Serotype),
        linetype = type,
        color = Serotype
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Symptomatic attack rate (%)") +
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(
      labels  = c(
        "Placebo",
        "Vaccine"
      )) +
    facet_wrap(~ Serostatus)
  
# plot symp attack rate by age and trial arm -----------------------
  
AR_plot_VJ = AR_model %>%
  filter(group == "AR_VJ") %>%
  separate(name, into = c("Arm", "Age")) %>% 
  mutate(Arm = factor(Arm, levels = 1:2, labels = c("placebo", "vaccine")),
          Age = factor(Age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"),
                        levels = 1:3)) %>% 
    bind_rows(AR_age_data) %>%
    ggplot(aes(x = Arm, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = Age,
        group = interaction(type, Age)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, Age),
        linetype = type,
        color = Age
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Symptomatic attack rate (%)") +
    scale_color_brewer(palette = "Accent") +
    guides(shape = "none",
           linetype = "none")


ggsave(
  plot = cowplot::plot_grid(
    AR_plot_BVK,
    AR_plot_VJ,
    labels = c("a", "b"),
    ncol = 1),
  filename = paste0(file_path, "/AR.png"),
  height = 60,
  width = 70,
  units = "cm",
  dpi = 300,
  scale = 0.75)


}

# function to plot vaccine efficacy --------------------------------------------

plot_BUT_VE  = function(file_path,
                    VE = NULL,
                    include_beta) {
  
VE_model = extract_BUT_model_results(VE)  
  
# plot VE by serostatus for each serotype --------------------------------------

VE_BKJT =  VE_model %>%
  separate(name, into = c("Serostatus","Serotype","Age", "Month")) %>% 
  mutate(Serostatus = factor(Serostatus, levels = 1:3,
                             labels = c("seronegative", "monotypic", "multitypic")),
         Age = factor(Age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"),
                      levels = 1:3),
         Serotype = factor(Serotype, levels = 1:2, labels = c(paste0("DENV", 1:2)))) %>% 
    mutate(Month = as.numeric(Month)) 
  
VE_plot = VE_BKJT %>% 
  filter(Serostatus != "multitypic") %>% 
  ggplot(aes(x = Month , y = mean)) +
  geom_line(aes(color = Serotype)) +
  geom_ribbon(aes( ymin = lower, ymax = upper,
        fill = Serotype), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine Efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 24, 12)) +
  facet_grid(Serostatus~Age) + 
  theme_light() +
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") 
  
  ggsave(
    plot = VE_plot,
    filename = paste0(file_path,  "/VE.png"),
    height = 21,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.8
  )
  
if(include_beta == 0){
  VE_plot2 = VE_BKJT %>% 
    filter(Serostatus != "multitypic") %>% 
    filter(Age == "2-6yrs") %>%  
    ggplot(aes(x = Month , y = mean)) +
    geom_line(aes(color = Serotype)) +
    geom_ribbon(aes( ymin = lower, ymax = upper,
                     fill = Serotype), alpha = 0.5) +
    labs(x = "Month", y = "Vaccine Efficacy (%)") +
    scale_x_continuous(breaks = seq(0, 24, 12)) +
    facet_grid(~Serostatus) + 
    theme_light() +
    theme(legend.position = "top")+
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired") 
  
  ggsave(
    plot = VE_plot2,
    filename = paste0(file_path,  "/VE_no_age.png"),
    height = 21,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.8
  )  
 
}
  
}


# function to plot titres ----------------------------------------------------------------------------------------
plot_BUT_titres = function(file_path, n = NULL){
  
  
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
           serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
           serotype = paste0("DENV", serotype))
  
  
plot_titres = model_titres %>%
  ggplot(aes(x = time, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper,), alpha = 0.2) +
  theme(legend.position = "top") +
  facet_grid(serostatus ~ serotype) +
  labs(x = "Month", y = "log titre")
  
  
ggsave(
  plot_titres,
  filename = paste0(file_path, "/plot_titres.jpg"),
  width = 40,
  height = 20,
  unit = "cm"
)

}

# function to plot everything --------------------------------------------------
plot_BUT_output = function(age_cases,
                           serotype_serostatus_cases,
                           file_path,
                           include_beta,
                           AR = NULL,
                           VE = NULL,
                           n = NULL) {
  
  library(tidyverse)
  library(Hmisc)
  library(cowplot)

  theme_set(
    theme_light() +
      theme(
        text = element_text(size = 16),
        legend.position ="top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)
      ))
  
  if (is.null(AR)) {
    AR = readRDS(paste0(file_path, "/AR.RDS"))
  }
  if (is.null(VE)) {
    VE = readRDS(paste0(file_path, "/VE.RDS"))
  }
  if (is.null(n)) {
    n = readRDS(paste0(file_path, "/n.RDS"))
  }

# plot attack rates 
plot_BUT_attack_rate(
  age_cases = age_cases,
  serotype_serostatus_cases = serotype_serostatus_cases,
  file_path = file_path,
  AR = AR,
  )
  
# plot VE
plot_BUT_VE(file_path = file_path, 
        VE = VE, 
        include_beta=include_beta)
  
# plot VE
plot_BUT_titres(file_path = file_path, n = n)

}
