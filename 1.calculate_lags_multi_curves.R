

# Figure 1: simulated data without any noise
# general params
source("~/Library/CloudStorage/OneDrive-UniwersytetJagielloński/Lags_part2/code/lag_calc_methods_comparison/0.config.R")
time.interval = 0.1
times = seq(0,Max.Time,time.interval)
lag = 2.5
growth.rate = 0.1
K = 5*10^6


# logistic model params
growth.rate.logistic = 0.25

# Monod model params: those parameters give as a curve similar to the logistic one
a = 0.03*10^7 
Vh = 150*10^(-7) 
Kh = 300
N0 = 10^6

# simulate data
byranayi_and_roberts.simulated.data.raw = data.frame(time = times,
                                                     biomass = Baranyi.Solution(
                                                       times, 
                                                       LOG10Nmax =  log10(K), 
                                                       mumax =  growth.rate, 
                                                       LOG10N0 = log10(N0), 
                                                       lag = lag))
logistic.simulated.data = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, lag, times)
monod.simulated.data = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=lag, times=times)
exponential.growth.simulated.data = Simulate.Exponential.Data.With.Lag(N0, growth.rate, lag, times)


all.simulated.data =
  rbind(byranayi_and_roberts.simulated.data.raw %>% mutate(curve_id = "Baranyi&Roberts"),
        logistic.simulated.data %>% mutate(curve_id = "Logistic"),
        monod.simulated.data %>% mutate(curve_id = "Monod"),
        exponential.growth.simulated.data %>% select(time, biomass) %>% mutate(curve_id = "Exponential"))

data.all.with.lag = suppressMessages(get_all_methods_lag(all.simulated.data, biomass.increase.threshold)) %>%
  rename(model = curve_id, method = lag_calculation_method)

data.all.with.lag$method <- factor(data.all.with.lag$method, levels = method.labels)



fig1 =
  ggplot(data.all.with.lag, aes(x = time, y = biomass)) +
  geom_rect(aes(xmin = -Inf, xmax = lag, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.1) +
  geom_line() +
  geom_line(aes(x = lag), color = "darkgrey") +
  scale_y_continuous(labels = label_scientific()) +
  facet_grid(cols = vars(method), 
             rows = vars(model),
             scales = "free",
             labeller = labeller(method = c("biomass increase" = "Biomass Increase",
                                            "max growth acceleration" = "Max Growth\nAcceleration",
                                            "tangent to max growth point" = "Tangent\nto the point",
                                            "tangent to max growth line" = "Tangent\nto the line",
                                            "par. fitting to baranyi model" = "Parameter Fitting\nto the Baranyi model",
                                            "par. fitting to logistic model" = "Parameter Fitting\nto the logistic model"))) +
  labs(x = "Time [h]", y = "CFU/mL") +
  theme_bw()+
  geom_text(aes(x = Inf, y = Inf, label = paste0("LAG = ", lag)), hjust = 2, vjust = 2, size = 3)+
  coord_cartesian(xlim = c(1, 23))

ggsave(fig1, 
       filename = sprintf("%sfig1.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 150, dpi = 300)




# Figure 2: experimental curves
fresh = read.delim(sprintf("%sfresh_culture.txt", EXPERIMENTAL.DATA.PATH)) %>% 
  mutate(curve_id = "Standard Curve")
colnames(fresh) = c("time", "biomass", "curve_id")

atypical = read.delim(sprintf("%schosen_exampless_biomass_ml.txt", EXPERIMENTAL.DATA.PATH)) %>% 
  filter(curve_id %in% c("curve_3", "curve_6", "curve_11")) %>%
  select(time = time_hours, biomass, curve_id) %>%
  mutate(type = case_when(
    curve_id == "curve_3" ~ "Atypical curve 1",
    curve_id == "curve_6" ~ "Atypical curve 2",
    curve_id == "curve_11" ~ "Atypical curve 3"
  )) %>%
  select(-curve_id) %>%
  rename(curve_id = type)

no.lag = read.delim(sprintf("%sBY_a_lag_data.txt", EXPERIMENTAL.DATA.PATH)) %>%
  filter(well == "E3") %>%
  select(time = time_h, biomass) %>%
  mutate(curve_id = "No-Lag")

all.data.fig2 = rbind(fresh, atypical, no.lag)
all.data.fig2.with.lag = suppressMessages(get_all_methods_lag(all.data.fig2, biomass.increase.threshold)) %>%
  rename(method = lag_calculation_method, type = curve_id)

all.data.fig2.with.lag$method <- factor(all.data.fig2.with.lag$method, levels = method.labels)


all.data.fig2.with.lag$type <- factor(all.data.fig2.with.lag$type, levels = c(
  "Standard Curve",
  "Atypical curve 1",
  "Atypical curve 2",
  "Atypical curve 3",
  "No-Lag"
))

fig2 = 
  ggplot(all.data.fig2.with.lag, aes(x = time, y = biomass)) +
  geom_rect(aes(xmin = -Inf, xmax = lag, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.9) +
  geom_line() +
  geom_line(aes(x = lag), color = "darkgrey") +
  scale_y_continuous(labels = label_scientific()) +
  facet_grid(cols = vars(method), 
             rows = vars(type),
             scales = "free_y",
             labeller = labeller(method = c("biomass increase" = "Biomass Increase",
                                            "max growth acceleration" = "Max Growth\nAcceleration",
                                            "tangent to max growth point" = "Tangent\nto the point",
                                            "tangent to max growth line" = "Tangent\nto the line",
                                            "par. fitting to baranyi model" = "Parameter Fitting\nto the Baranyi model",
                                            "par. fitting to logistic model" = "Parameter Fitting\nto the logistic model"))) +
  labs(x = "Time", y = "CFU/mL") +
  theme_bw()+
  geom_text(aes(x = Inf, y = Inf, label = paste0("LAG = ", lag)), hjust = 1.5, vjust = 6, size = 3)+
  coord_cartesian(xlim = c(1, 23))

ggsave(fig2, 
       filename = sprintf("%sfig2.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 250, height = 180, dpi = 300)






# Supplementary Figures
# S1 Fitting models to nolag data
nolag.data.av = no.lag %>% 
  group_by(time) %>% 
  summarise(biomass = mean(biomass)) %>% 
  ungroup() %>% 
  mutate(curve_id = "no_lag_curve") %>% 
  arrange(time) 

# Fit no lag data to models
pars = get_def_pars()
pars$model <- "baranyi"
calc_lag(nolag.data.av, method = "parameter fitting to a model", pars) %>% distinct(curve_id, lag)
nlsmod = Fit.To.Baranyi(nolag.data.av,
                        max_iter = 100, 
                        init_lag = NULL, 
                        init_gr_rate = NULL)

# Note here we get a different lag than using lag calculator
nlsmod.log = Fit.To.Logistic.With.Lag(N0=nolag.data.av$biomass[1],
                                      data = nolag.data.av %>% select(time, biomass),
                                      lower.bound = c(0.01,1,0), 
                                      upper.bound = c(10,10^10,10))


log.pars.fitted = nlsmod.log$optim$bestmem
logistic.model.params = data.frame(growth.rate = log.pars.fitted[1] %>% as.numeric(),
                                   K = log.pars.fitted[2] %>% as.numeric(),
                                   lag = log.pars.fitted[3] %>% as.numeric(),
                                   N0 = nolag.data.av$biomass[1])
baranyi.model.params = data.frame(
  LOG10Nmax =  coef(nlsmod)[names(coef(nlsmod)) == "LOG10Nmax"] %>% as.numeric(), 
  mumax =  coef(nlsmod)[names(coef(nlsmod)) == "mumax"] %>% as.numeric(), 
  LOG10N0 = coef(nlsmod)[names(coef(nlsmod)) == "LOG10N0"] %>% as.numeric(), 
  lag = coef(nlsmod)[names(coef(nlsmod)) == "lag"] %>% as.numeric())

nolag.data.av$logistic = Simulate.Logistic.With.Lag(N0=nolag.data.av$biomass[1],
                                                    growth.rate = logistic.model.params$growth.rate,
                                                    K = logistic.model.params$K,
                                                    lag = logistic.model.params$lag,
                                                    times = nolag.data.av$time) %>%
  pull(biomass)


nolag.data.av$Baranyi = Baranyi.Solution(nolag.data.av$time, 
                                         LOG10Nmax = baranyi.model.params$LOG10Nmax, 
                                         mumax =  baranyi.model.params$mumax,
                                         LOG10N0 = baranyi.model.params$LOG10N0, 
                                         lag = baranyi.model.params$lag)


g=ggplot(nolag.data.av %>% 
           rename(`real CFU` = biomass) %>%
           tidyr::gather(key = "Model", value = "CFU", 
                         `real CFU`, Baranyi, logistic) %>%
           mutate(`log10(CFU)`= log10(CFU)) %>%
           tidyr::gather(key = "scale", value = "Cell.count",CFU, `log10(CFU)`),
         aes(x=time, y = Cell.count, col = Model, linetype = Model, alpha = Model)) + 
  geom_line(size = 1) +
  theme_bw() +
  facet_grid(scale ~ ., scales = "free") +
  scale_linetype_manual(values = c("real CFU" = "solid", "Baranyi" = "dashed", logistic = "dashed")) +
  scale_color_manual(values = c("real CFU" = "darkgray", "Baranyi" = "darkgreen", logistic = "darkviolet")) +
  scale_alpha_manual(values = c("real CFU" = 1, "Baranyi" = 0.75, logistic = 0.75)) +
  xlab("Time [h]") + 
  ylab("")
g
ggsave(sprintf("%sS1.png", OUTPUT_FIGS_PATH), g, width = 16, height=8, units = "cm")



# THESE WILL BE USED TO GENERATE FIG 3
saveRDS(logistic.model.params, 
        paste0(OUTPUT_FIGS_PATH,"logistic.model.params.rds"))
saveRDS(baranyi.model.params, 
        paste0(OUTPUT_FIGS_PATH,"baranyi.model.params.rds"))



# S2a
curves11 = read.delim(sprintf("%schosen_exampless_biomass_ml.txt", EXPERIMENTAL.DATA.PATH))
curve_names = c("curve_1",
                "curve_2",
                "curve_3",
                "curve_4",
                "curve_5",
                "curve_6",
                "curve_7",
                "curve_8",  
                "curve_9",
                "curve_10",
                "curve_11")
curves11$curve_id = factor(curves11$curve_id,
                           levels = curve_names)

supfig2 = 
  ggplot(curves11, aes(x = time_hours, y = biomass), size = 0.5) +
  geom_point() +
  geom_line() +
  facet_wrap("curve_id",
             nrow = 6,
             ncol = 2,
             scales = "free_y") +
  theme_bw() +
  xlab("Time [h]") +
  ylab("CFU") +
  ggtitle("A")

ggsave(supfig2,
       filename = sprintf("%ssup.fig2a.jpg", OUTPUT_FIGS_PATH),
       units = "mm", width = 200, height = 150, dpi = 300)


# S2b
real.data = curves11 %>%
  select(time = time_hours, biomass = biomass_ml, curve_id)
real.data.with.lag = suppressMessages(get_all_methods_lag(real.data, biomass.increase.threshold))

data.smooth = smooth_data(real.data)
real.data.smooth.with.lag =  suppressMessages(get_all_methods_lag(data.smooth, biomass.increase.threshold))

data.short = cut_the_data(real.data,12)
real.data.short.with.lag = suppressMessages(get_all_methods_lag(data.short, biomass.increase.threshold))

lag.data =
  rbind(real.data.smooth.with.lag %>%
          distinct(lag, curve_id, lag_calculation_method) %>%
          mutate(data.type = "smoothened"),
        real.data.with.lag %>%
          distinct(lag, curve_id, lag_calculation_method) %>%
          mutate(data.type = "original")) %>%
  rbind(
    real.data.short.with.lag %>%
      distinct(lag, curve_id, lag_calculation_method) %>%
      mutate(data.type = "cut")) %>%
  mutate(curve_id = factor(curve_id, levels = curve_names)) %>%
  mutate(Method = plyr::mapvalues(lag_calculation_method,
    from = c("biomass increase",  "max growth acceleration",  "tangent to max growth point",
             "tangent to max growth line", "par. fitting to baranyi model",  "par. fitting to logistic model"),
    to = c("Biomass Increase",  "Max Growth\nAcceleration", "Tangent\nto the point","Tangent\nto the line",
           "Parameter Fitting\nto the Baranyi model","Parameter Fitting\nto the logistic model")))

g=ggplot(lag.data %>%
           rename(`Data type` = data.type),
         aes(y = lag, x = curve_id)) +
  geom_jitter(aes(col = Method), height = 0, size = 2.5) +
  geom_boxplot(alpha = 0.2, outliers = FALSE) +
  facet_grid(. ~`Data type`) +
  Get.Theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom") +
  xlab("") +
  ylab("Estimated lag") +
  ggtitle("B")
g
ggsave(sprintf("%ssup.fig.2b.png", OUTPUT_FIGS_PATH), g, width = 200, height=120, units = "mm")

