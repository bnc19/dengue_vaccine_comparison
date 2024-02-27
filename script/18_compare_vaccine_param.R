# Create a figure comparing the parameters of the best fitting models for each vaccine
library(tidyverse)

model_dep = data.frame(
vaccine = c("Dengvaxia", "Qdenga", "Butantan-DV"), 
L = c("serostatus & serotype", "serostatus", NA),
delta = c("age", "serotype", NA),
lc = c("serostatus & serotype (serostatus = MO, MU)", 
       "serostatus & serotype (serostatus = MO, MU)",
       "serostatus & serotype (serostatus = MO, MU)"), 
beta = c(NA, "age (youngest)", NA),
p = c("age", "age", "age"),
tau = c("global", "serotype", NA ), 
omega = c("global", "global", "global"),
kappa = c(NA, NA, NA),
eta = c(NA, NA, "age (youngest)"),
rho = c("global", "global", "global"),
w = c("global", "global", "global"),
alpha = c("global", "global", NA)
)

# format data ready to plot 
df = model_dep %>%
  pivot_longer(cols = L:alpha,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "serotype",
        "age",
        "age (youngest)",
        "serostatus",
        "serostatus & serotype (serostatus = MO, MU)",
        "serostatus & serotype",
        "global",
        "not included"
      ) ))

mycols = c("#A6BDDB", "#1C9099", "#6BAED6",
           "#666699", "#FBB4B9", 
           "#FEEBE2", "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1 = df %>%
  ggplot(aes(x = Parameter, y = vaccine)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_fill_manual(values=mycols) +
  scale_y_discrete(labels = c('omega' = expression(omega),
                              "kappa" = expression(kappa),
                              "beta" = expression(beta),
                              "delta" = expression(delta),
                              "alpha" = expression(alpha),
                              "tau" = expression(tau),
                              "rho" = expression(rho),
                              "eta" = expression(eta))) +
  ylab("Vaccine")




ggsave(p1, file = "compare_vaccines/output/compare_model_param.jpg",
       height = 20, width = 30, units="cm", scale = 0.7)

