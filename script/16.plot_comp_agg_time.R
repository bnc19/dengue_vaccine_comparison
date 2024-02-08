# Plot 2 year VE for all vaccines by serostatus and serotype (but not age)
rm(list = ls())

source("compare_vaccines/R/plotting_functions.R")
library(tidyverse)
n= 100

TAK_file = "TAK/output/M2/VE.RDS"
CYD_file = "CYD/output/M12/VE.RDS"
BUT_file = "BUT/output/M4/VE.RDS"

# Read in best fitting models and sample n iterations 

BUT_VE = as.data.frame(sapply(readRDS(BUT_file), sample, n))

# Select VE without age 
CYD_VE = select(readRDS(CYD_file),contains("VE_BKRT"))
TAK_VE = select(readRDS(TAK_file),contains("VE_BKRT"))

CYD_VE2 = as.data.frame(sapply(CYD_VE, sample, n))
TAK_VE2 = as.data.frame(sapply(TAK_VE, sample, n))

# summarise mean and 95% CrI VE by serostatus, serotype at 24 months
tidy_BUT_VE = format_2_year_VE(fit=BUT_VE, vaccine = "BUT")
tidy_CYD_VE = format_2_year_VE(fit=CYD_VE2, vaccine = "CYD")
tidy_TAK_VE = format_2_year_VE(fit=TAK_VE2, vaccine = "TAK")


tidy_BUT_VE$Vaccine = "Butantan-DV"
tidy_CYD_VE$Vaccine = "Dengvaxia"
tidy_TAK_VE$Vaccine = "Qdenga"

year2_comp = tidy_TAK_VE %>% 
  bind_rows(tidy_BUT_VE, tidy_CYD_VE) %>%  
  ggplot(aes(x = Vaccine, y= mean, group = Serotype)) +
  geom_point(aes(color = Serotype),
             position = position_dodge(0.4), size = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Serotype),
                position = position_dodge(0.4), width = 0.5) +
  facet_wrap(~Serostatus, ncol =3) + 
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    legend.position = c(0.94,0.25),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(0.2, "cm")
  )  + ylab("Vaccine efficacy (%)") +
  ggtitle("Cumulative vaccine efficacy across months 1 to 24")

ggsave(
  plot = year2_comp,
  filename = "compare_vaccines/output/year_2_VE_comp.jpg",
  height = 10,
  width = 30,
  units = "cm",
  dpi = 600,
  scale = 0.6
)

# summarise mean and 95% CrI VE by outcome, seerostatus, serotype upto 4.5 yrs
tidy_CYD_VE_4.5 = format_4.5_year_VE(CYD_VE2)
tidy_TAK_VE_4.5 = format_4.5_year_VE(TAK_VE2)

tidy_CYD_VE_4.5$Vaccine = "Dengvaxia"
tidy_TAK_VE_4.5$Vaccine = "Qdenga"

year_4.5_comp = tidy_TAK_VE_4.5 %>% 
  bind_rows(tidy_CYD_VE_4.5) %>%  
  ggplot(aes(x = Vaccine, y= mean, group = Serotype)) +
  geom_point(aes(color = Serotype),
             position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Serotype),
                position = position_dodge(0.4), width = 0.5) +
  facet_grid(Outcome~Serostatus) + 
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    legend.position ="top",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )  + ylab("Vaccine efficacy (%)") +
  ggtitle("Cumulative vaccine efficacy across months 1 to 54")

ggsave(
  plot = year_4.5_comp,
  filename = "output/year_4.5_comp.jpg",
  height = 20,
  width = 30,
  units = "cm",
  dpi = 600,
  scale = 0.6
)

