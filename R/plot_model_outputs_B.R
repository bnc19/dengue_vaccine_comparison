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

plot_BUT_attack_rate = function(cases,
                                file_path,
                                AR) {
  
age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill = c(scales::brewer_pal(palette = "RdPu")(6)[2:5]) 
  
# add aggregated populations to data and calculate attack rates
AR_model = extract_BUT_model_results(AR)

BVK_cases = calc_BUT_attack_rates(cases$Sy_BVK)
  
# plot serotype serostatus attack rate 
AR_plot_BVK = AR_model %>%
  filter(group == "AR_BVK") %>%
  separate(name, into = c("serostatus", "arm", "serotype")) %>% 
  mutate(arm = factor(arm, levels = 1:2, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, levels = 1:2, labels = c("seronegative", "seropositive")),
         serotype = factor(serotype, levels = 1:2, labels = c(paste0("DENV", 1:2)))) %>% 
    bind_rows(BVK_cases) %>%
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
    labs(x = " ", y = "Symptomatic attack rate (%)") +
    scale_color_manual(values = serotype_fill) +
    facet_wrap( ~serostatus)
  
# plot symp attack rate by age and trial arm -----------------------
  
BVJ_cases = calc_BUT_attack_rates(cases$Sy_BVJ)

AR_plot_BVJ = AR_model %>%
  filter(group == "AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(arm = factor(arm, levels = 1:2, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, levels = 1:2, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"),
                        levels = 1:3)) %>% 
    bind_rows(BVJ_cases) %>%
    ggplot(aes(x = arm, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = age,
        group = interaction(type, age)
      ),
      position = position_dodge(width = 0.5),
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
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "Symptomatic attack rate (%)") +
    scale_color_brewer(palette = "Accent") +
    guides(shape = "none",
           linetype = "none") +
  facet_wrap(~serostatus) +
  scale_color_manual(values = age_fill) 
  


ggsave(
  plot = cowplot::plot_grid(
    AR_plot_BVK,
    AR_plot_BVJ,
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
  separate(name, into = c("serostatus","serotype","age", "month")) %>% 
  mutate(serostatus = factor(serostatus, levels = 1:3,
                             labels = c("seronegative", "monotypic", "multitypic")),
         age = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"),
                      levels = 1:3),
         serotype = factor(serotype, levels = 1:2, labels = c(paste0("DENV", 1:2)))) %>% 
    mutate(month = as.numeric(month)) 
  
VE_plot = VE_BKJT %>% 
  filter(serostatus != "multitypic") %>% 
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = serotype)) +
  geom_ribbon(aes( ymin = lower, ymax = upper,
        fill = serotype), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine Efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 24, 12)) +
  facet_grid(serostatus~age) + 
  theme_light() +
  theme(legend.position = "top") + 
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
    filter(serostatus != "multitypic") %>% 
    filter(age == "2-6yrs") %>%  
    ggplot(aes(x = Month , y = mean)) +
    geom_line(aes(color = serotype)) +
    geom_ribbon(aes( ymin = lower, ymax = upper,
                     fill = serotype), alpha = 0.5) +
    labs(x = "Month", y = "Vaccine efficacy (%)") +
    scale_x_continuous(breaks = seq(0, 24, 12)) +
    facet_grid(~serostatus) + 
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
plot_BUT_output = function(cases,
                           file_path,
                           include_beta) {
  
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
  

    AR = readRDS(paste0(file_path, "/AR.RDS"))
    VE = readRDS(paste0(file_path, "/VE.RDS"))
    n = readRDS(paste0(file_path, "/n.RDS"))

# plot attack rates 
plot_BUT_attack_rate(
  cases = cases,
  file_path = file_path,
  AR = AR
  )
  
# plot VE
plot_BUT_VE(file_path = file_path,
            VE = VE,
            include_beta = include_beta)
  
# plot VE
plot_BUT_titres(file_path = file_path, n = n)

}
