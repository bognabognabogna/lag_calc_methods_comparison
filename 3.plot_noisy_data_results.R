methods_labeller = c("biomass increase" = "Biomass Increase",
                     "max growth acceleration" = "Max Growth\nAcceleration",
                     "tangent to max growth point" = "Tangent\nto the point",
                     "tangent to max growth line" = "Tangent\nto the line",
                     "par. fitting to baranyi model" = "Parameter Fitting\nto the Baranyi model",
                     "par. fitting to logistic model" = "Parameter Fitting\nto the logistic model")
sd_range = c(0.01, 0.1, 1)

# Figure 3 PLOT
## GROWTH
growth = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.growth.rate.rds")) %>%
  mutate(real.minus.estimated.lag = real.lag - lag)
growth$sd = as.factor(growth$sd)
growth$lag_calculation_method = as.factor(growth$lag_calculation_method)
growth$lag_calculation_method = ordered(growth$lag_calculation_method, levels = method.labels)
growth <- growth %>% mutate(GR = paste0("GR = ", round(growth.rate, 2)))

growth.stats = growth %>%
  group_by(sd, lag_calculation_method, GR, time.interval, real.lag) %>%
  summarise(
    mean.lag = mean(lag, na.rm = TRUE),
    n = n_distinct(curve_id),
    `bias [h]` = mean.lag - unique(real.lag),
    # variance i.e. precision
    `var [h^2]` = mean((lag - mean.lag)^2, na.rm = TRUE)) %>%
  tidyr::gather(key = "Measure", value = "value",  `bias [h]`, `var [h^2]`)
    

plot_growth_stats = 
  ggplot(growth.stats %>% 
           ungroup() %>%
           mutate(sd = as.numeric(as.character(sd))),
         aes(x = sd, y = value, col = GR))+
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dashed") +
  facet_grid(cols = vars(lag_calculation_method),
             rows = vars(Measure),
             labeller = labeller(lag_calculation_method = methods_labeller),
             scales = "free") +
  theme_bw()+
  scale_x_continuous(breaks=sd_range, labels = sd_range, name = "Noise [sd]", trans = 'log10') +
  ylab(" ")+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 1, legend.position = "bottom")+
  scale_color_manual(values = c("GR = 0.71" = "lightblue", "GR = 1.42" = "blue", "GR = 2.83" = "darkblue"), name = "Growth rate")

plot_growth_stats
ggsave(plot_growth_stats,
       filename = sprintf("%sfig3b_stats.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 140, dpi = 300)






plot_growth = 
  ggplot(growth, #%>% #filter(lag_calculation_method != "par. fitting to baranyi model"), 
         aes(x = sd, y = obs.minus.real.lag, col = lag_calculation_method))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.25, 
               aes(fill = lag_calculation_method),
               outlier.shape = NA)+
  geom_jitter(
    width = 0.1,
    alpha = 0.3)+
  facet_grid(cols = vars(lag_calculation_method),
             rows = vars(GR),
             labeller = labeller(lag_calculation_method = methods_labeller),
             scales = "free_y") +
  theme_bw()+
  stat_summary(geom = "point",
               fun = "median",
               size = 2,
               shape = 15,
               colour = "grey")+
  xlab("Noise [sd]")+
  ylab("Estimated - real lag")+
  theme(legend.position = "none")+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("#F59A22", "#BC4198", "#9FBBE5", "#2C5AA0", "#63E383", "#1C9C3C"))+
  scale_color_manual(values = c("#F59A22", "#BC4198", "#9FBBE5", "#2C5AA0", "#63E383", "#1C9C3C"))+
  scale_y_continuous(limits = c(1.1*min(growth$obs.minus.real.lag, na.rm = T), 10)) +
  ggtitle("B. Growth Rate")

plot_growth

ggsave(plot_growth,
       filename = sprintf("%sfig3b.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 140, dpi = 300)


# TIME INTERVAL

interval = 
  readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.time.interval.rds"))
interval$sd = as.factor(interval$sd)
interval$lag_calculation_method = as.factor(interval$lag_calculation_method)
levels(interval$lag_calculation_method)
interval$lag_calculation_method = ordered(interval$lag_calculation_method, levels = method.labels)

interval.stats = interval %>%
  group_by(sd, lag_calculation_method, growth.rate, time.interval, real.lag) %>%
  summarise(
    mean.lag = mean(lag, na.rm = TRUE),
    n = n_distinct(curve_id),
    `bias [h]` = mean.lag - unique(real.lag),
    # vaiance i.e. precision
    `var [h^2]` = mean((lag - mean.lag)^2, na.rm = TRUE)) %>%
  tidyr::gather(key = "Measure", value = "value", `bias [h]`, `var [h^2]`)



plot_interval_stats = 
  ggplot(interval.stats %>%
           mutate(sd = as.numeric(as.character(sd)),
                  time.interval = paste0(time.interval, "h"), 
                  time.interval = as.factor(time.interval)) %>%
                    # HACK ALERT! MENTION THAT IN MANUSCRIPT
            mutate(value = if_else(value > 15, NA, value)),
         aes(x = sd, y = value, col = time.interval))+
  geom_hline(yintercept = 0, cool = "black", linetype = "dashed")+
  geom_point() +
  geom_line() +
  facet_grid(cols = vars(lag_calculation_method),
             rows = vars(Measure),
             labeller = labeller(lag_calculation_method = methods_labeller),
             scales = "free") +
  theme_bw()+
  scale_x_continuous(breaks=sd_range, labels = sd_range, name = "Noise [sd]", trans = 'log10') +
  ylab(" ")+
  #theme(legend.position = "none")+
  #coord_cartesian(ylim = c(0, 1.5)) +
  theme(aspect.ratio = 1, legend.position = "bottom")+
  scale_color_manual(values = c("0.1h" = "darkblue", "0.5h" = "blue", "1h" = "lightblue"), name = "Measurement interval")


plot_interval_stats

ggsave(plot_interval_stats,
       filename = sprintf("%sfig3a_stats.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 140, dpi = 300)




plot_interval = 
  ggplot(interval, #%>% filter(lag_calculation_method != "par. fitting to baranyi model"), 
         aes(x = sd, y = obs.minus.real.lag, col = lag_calculation_method))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.25, 
               aes(fill = lag_calculation_method),
               outlier.shape = NA)+
  geom_jitter(
    width = 0.1,
    alpha = 0.3)+
  
  facet_grid(cols = vars(lag_calculation_method),
             rows = vars(time.interval),
             labeller = labeller(time.interval = c("0.1" = "TI 0.1h", 
                                                   "0.5" = "TI 0.5h", 
                                                   "1" = "TI 1h"),
                                 lag_calculation_method = methods_labeller),
             scales = "free_y")+
  theme_bw()+
  stat_summary(geom = "point",
               fun = "median",
               size = 2,
               shape = 15,
               colour = "grey")+
  xlab("Noise [sd]")+
  ylab("Estimated - real lag")+
  theme(legend.position = "none")+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("#F59A22", "#BC4198", "#9FBBE5", "#2C5AA0", "#63E383", "#1C9C3C"))+
  scale_color_manual(values = c("#F59A22", "#BC4198", "#9FBBE5", "#2C5AA0", "#63E383", "#1C9C3C")) +
  ggtitle("A. Measurement frequency")

plot_interval

ggsave(plot_interval,
       filename = sprintf("%sfig3a.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 140, dpi = 300)




## real lag
lagvary = 
  readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.lags.rds")) %>%
  mutate(RL = round(real.lag, 2))

lagvary$sd = as.factor(lagvary$sd)
lagvary$lag_calculation_method = as.factor(lagvary$lag_calculation_method)
levels(lagvary$lag_calculation_method)
lagvary$lag_calculation_method = ordered(lagvary$lag_calculation_method, levels = method.labels)


lagvary.stats = lagvary %>%
  group_by(sd, lag_calculation_method, growth.rate, time.interval, RL) %>%
  summarise(
    mean.lag = mean(lag, na.rm = TRUE),
    n = n_distinct(curve_id),
    `bias [h]` = mean.lag - unique(real.lag),
    # vaiance i.e. precision
    `var [h^2]` = mean((lag - mean.lag)^2, na.rm = TRUE)) %>%
  tidyr::gather(key = "Measure", value = "value", `bias [h]`, `var [h^2]`)


plot_lagvary_stats = 
  ggplot(lagvary.stats %>%
           mutate(sd = as.numeric(as.character(sd)),
                  RL = as.character(RL)),
         aes(x = sd, y = value, col = RL))+
  geom_hline(yintercept = 0, col = "black", linetype = "dashed")+
  geom_point() +
  geom_line() +
  facet_grid(cols = vars(lag_calculation_method),
             rows = vars(Measure),
             labeller = labeller(lag_calculation_method = methods_labeller),
             scales = "free") +
  theme_bw()+
  scale_x_continuous(breaks=sd_range, labels = sd_range, name = "Noise [sd]", trans = 'log10') +
  ylab(" ")+
  theme(aspect.ratio = 1, legend.position = "bottom")+
  scale_color_manual(values = c("0.25" = "darkblue", "1" = "blue", "2.5" = "lightblue"), name = "Real Lag Duration")
plot_lagvary_stats

ggsave(plot_lagvary_stats,
       filename = sprintf("%sfig3c_stats.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 140, dpi = 300)




plot_lagvary = 
  ggplot(lagvary, #%>% filter(lag_calculation_method != "par. fitting to baranyi model"), 
         aes(x = sd, y = obs.minus.real.lag, col = lag_calculation_method))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.25, 
               aes(fill = lag_calculation_method),
               outlier.shape = NA)+
  geom_jitter(
    width = 0.1,
    alpha = 0.3)+
  facet_grid(cols = vars(lag_calculation_method),
             rows = vars(RL),
             labeller = labeller(
               RL = c("0.009" = "Lag 0.009h", 
                      "1" = "Lag 1h", 
                      "2.5" = "Lag 2.5h"),
               lag_calculation_method = methods_labeller)) +
  theme_bw()+
  stat_summary(geom = "point",
               fun = "median",
               size = 2,
               shape = 15,
               colour = "grey")+
  xlab("Noise [sd]")+
  ylab("Estimated - real lag")+
  theme(legend.position = "none")+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("#F59A22", "#BC4198", "#9FBBE5", "#2C5AA0", "#63E383", "#1C9C3C"))+
  scale_color_manual(values = c("#F59A22", "#BC4198", "#9FBBE5", "#2C5AA0", "#63E383", "#1C9C3C")) +
  ggtitle("Varying lag durations")
plot_lagvary
ggsave(plot_lagvary,
       filename = sprintf("%ssup.fig4.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 140, dpi = 300)


fig3 = plot_interval_stats / plot_growth_stats / plot_lagvary_stats
fig3 = fig3 +  plot_annotation(tag_levels = 'A') +  theme(plot.title = element_text(size = 20))
ggsave(filename = sprintf("%sfig3.jpg", OUTPUT_FIGS_PATH), fig3, width = 12, height = 13, dpi = 600)
