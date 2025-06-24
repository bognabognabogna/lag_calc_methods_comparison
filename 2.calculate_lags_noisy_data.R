source("~/Library/CloudStorage/OneDrive-UniwersytetJagielloński/Lags_part2/code/lag_calc_methods_comparison/0.config.R")

logistic.model.params = readRDS(paste0(OUTPUT_FIGS_PATH,"logistic.model.params.rds"))
baranyi.model.params = readRDS(paste0(OUTPUT_FIGS_PATH,"baranyi.model.params.rds"))

########################### simulated data with noise ################
time.intervals = c(0.1, 0.5, 1)
Num.obs = 500
sd_range = c(0.01, 0.1, 1) 
time.interval = 0.1
times = seq(0,Max.Time,time.interval)



# LOGISTIC PARAMS
real.lag.logistic = logistic.model.params$lag
growth.rate.logistic = logistic.model.params$growth.rate
K = logistic.model.params$K
N0 = logistic.model.params$N0
# in case simulated data go below 0 change them to MIN.BIOMASS
MIN.BIOMASS = 1



# example data with noise
logistic.simulated.data.example.1 = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, real.lag.logistic, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = sd_range[1]*N0)) %>%
  mutate(sd = sd_range[1], K = K, real.lag = real.lag.logistic, growth.rate = growth.rate.logistic)
logistic.simulated.data.example.2 = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, real.lag.logistic, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = sd_range[2]*N0))%>%
  mutate(sd = sd_range[2], K = K, real.lag = real.lag.logistic, growth.rate = growth.rate.logistic)
logistic.simulated.data.example.3 = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, real.lag.logistic, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = sd_range[3]*N0))%>%
  mutate(sd = sd_range[3], K = K, real.lag = real.lag.logistic, growth.rate = growth.rate.logistic)

logistic.simulated.data.example = rbind(logistic.simulated.data.example.1,
                                        rbind(logistic.simulated.data.example.2,
                                              logistic.simulated.data.example.3)) %>%
  rowwise() %>%
  mutate(sd = paste0("sd = ", sd),
         biomass = max(biomass, MIN.BIOMASS))

jpeg(sprintf("%sFig_example_noisy_curve_logistic.png", OUTPUT_FIGS_PATH), width = 20, height=10, units = "cm", res = 600)
ggplot(logistic.simulated.data.example %>% filter(time <= 4), aes(x = time, y = biomass)) +
  #geom_point() + 
  geom_line(size = 0.5) +
  facet_grid(.~sd) + 
  theme_bw() +
  scale_y_continuous(trans='log10', name = "CFU/mL") +
  xlab("Time [h]")
dev.off()


################## LOGISTIC MODEL SIMULATIONS ###########################
growth.rates.logistic= growth.rate.logistic*c(0.5, 1, 2)
real.lags.logistic = c(real.lag.logistic, 1, 2.5)


all.lag.data.logistic = data.frame(curve_id = character(0),
                          lag = numeric(0),
                          lag_calculation_method = character(0),
                          sd = numeric(0),
                          growth.rate = numeric(0),
                          real.lag = numeric(0),
                          time.interval = numeric(0))

logistic.curves = data.frame(time = numeric(0),
                             biomass = numeric(0),
                             curve_id = character(0),
                             sd = numeric(0),
                             growth.rate = numeric(0),
                             real.lag = numeric(0),
                             time.interval = numeric(0))

all.lag.data.logistic.2 = all.lag.data.logistic %>% filter(FALSE)
logistic.curves.2 = logistic.curves %>% filter(FALSE)
for (this.real.lag in real.lags.logistic) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, this.real.lag, times)
  lag.df.obj = suppressMessages(Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                                           sd_range = sd_range,
                                                                           biomass.increase.threshold,
                                                                           Num.obs = Num.obs,
                                                                           MIN.BIOMASS= MIN.BIOMASS))
  
  logistic.curves.2 = rbind(logistic.curves.2,
                            lag.df.obj$curves %>% mutate(growth.rate = growth.rate.logistic,
                                                         real.lag = this.real.lag,
                                                         time.interval = time.interval))
  all.lag.data.logistic.2 = rbind(all.lag.data.logistic.2,
                         lag.df.obj$lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                                      real.lag = this.real.lag,
                                                      time.interval = time.interval))
}


all.lag.data.logistic.3 = all.lag.data.logistic %>% filter(FALSE)
logistic.curves.3 = logistic.curves %>% filter(FALSE)
for (this.growth.rate in growth.rates.logistic) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, this.growth.rate, K, real.lag.logistic, times)
  lag.df.obj = suppressMessages(Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                                           sd_range = sd_range,
                                                                           biomass.increase.threshold,
                                                                           Num.obs = Num.obs, 
                                                                           MIN.BIOMASS= MIN.BIOMASS))
  logistic.curves.3 = rbind(logistic.curves.3,
                            lag.df.obj$curves %>% mutate(growth.rate = this.growth.rate,
                                                         real.lag = real.lag.logistic,
                                                         time.interval = time.interval))
  all.lag.data.logistic.3 = rbind(all.lag.data.logistic.3,
                         lag.df.obj$lag.df %>% mutate(growth.rate = this.growth.rate,
                                                      real.lag = real.lag.logistic,
                                                      time.interval = time.interval))
}


all.lag.data.logistic.4 = all.lag.data.logistic %>% filter(FALSE) %>% mutate(time.interval = numeric(0))
logistic.curves.4 = logistic.curves %>% filter(FALSE)
for (this.time.interval in time.intervals) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0,
                                                             growth.rate.logistic,
                                                             K,
                                                             real.lag.logistic,
                                                             seq(0,Max.Time,this.time.interval))
  lag.df.obj = suppressMessages(Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                                           sd_range = sd_range,
                                                                           biomass.increase.threshold,
                                                                           Num.obs = Num.obs,
                                                                           MIN.BIOMASS= MIN.BIOMASS))
  logistic.curves.4 = rbind(logistic.curves.4,
                            lag.df.obj$curves %>% mutate(growth.rate = growth.rate.logistic,
                                                         real.lag = real.lag.logistic,
                                                         time.interval = this.time.interval))
  all.lag.data.logistic.4 = rbind(all.lag.data.logistic.4, lag.df.obj$lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                                                      real.lag = real.lag.logistic,
                                                                      time.interval = this.time.interval))
}

all.lag.data.logistic.2 = all.lag.data.logistic.2 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.logistic.3 = all.lag.data.logistic.3 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.logistic.4 = all.lag.data.logistic.4 %>% mutate(obs.minus.real.lag = lag - real.lag)

saveRDS(all.lag.data.logistic.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.lags.rds"))
saveRDS(all.lag.data.logistic.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.growth.rate.rds"))
saveRDS(all.lag.data.logistic.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.time.interval.rds"))
saveRDS(logistic.curves.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.lags.curves.rds"))
saveRDS(logistic.curves.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.growth.rate.curves.rds"))
saveRDS(logistic.curves.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.time.interval.curves.rds"))








# BARANYI

# baranyi params
growth.rate.baranyi = baranyi.model.params$mumax
LOG10Nmax = baranyi.model.params$LOG10Nmax
LOG10N0 = baranyi.model.params$LOG10N0
real.lag.baranyi = baranyi.model.params$lag


growth.rates.baranyi = growth.rate.baranyi*c(0.5, 1,2)
real.lags.baranyi = c(0, real.lag.baranyi, 2.5)

baranyi.simulated.data.basic = baranyi.simulated.data.basic = data.frame(time = times,
                                                                         biomass = Baranyi.Solution(times, LOG10Nmax, growth.rate.baranyi, LOG10N0, real.lag.baranyi))

baranyi.simulated.data.example.1 = baranyi.simulated.data.basic %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = sd_range[1]*N0))%>%
  mutate(sd = sd_range[1],real.lag = real.lag.baranyi, growth.rate = growth.rate.baranyi)
baranyi.simulated.data.example.2 = baranyi.simulated.data.basic %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = sd_range[3]*N0))%>%
  mutate(sd = sd_range[3], real.lag = real.lag.baranyi, growth.rate = growth.rate.baranyi)
baranyi.simulated.data.example.3 = baranyi.simulated.data.basic %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = sd_range[4]*N0))%>%
  mutate(sd = sd_range[4],  real.lag = real.lag.baranyi, growth.rate = growth.rate.baranyi)
baranyi.simulated.data.example = rbind(baranyi.simulated.data.example.1,
                                       rbind(baranyi.simulated.data.example.2,
                                             baranyi.simulated.data.example.3)) %>%
  rowwise() %>%
  mutate(sd = paste0("sd = ", sd))

write.csv(baranyi.simulated.data.example %>%
            select(time, biomass, curve_idd = sd), 
          file = sprintf("%sbaranyi.ssimulated_data.with.noise.csv", OUTPUT_FIGS_PATH), row.names = FALSE)


jpeg(sprintf("%sFig_example_noisy_curve_baranyi.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(baranyi.simulated.data.example, aes(x = time, y = biomass)) +
  geom_point() + geom_line() +
  facet_grid(.~sd) + 
  theme_bw() +
  xlab("time [h]") +
  ylab("biomass [CFU/mL]") 
dev.off()











all.lag.data.baranyi = data.frame(curve_id = character(0),
                          lag = numeric(0),
                          lag_calculation_method = character(0),
                          sd = numeric(0),
                          growth.rate = numeric(0),
                          real.lag = numeric(0),
                          time.interval = numeric(0))
baranyi.curves = data.frame(time = numeric(0),
                             biomass = numeric(0),
                             curve_id = character(0),
                             sd = numeric(0),
                             growth.rate = numeric(0),
                             real.lag = numeric(0),
                             time.interval = numeric(0))

all.lag.data.baranyi.2 = all.lag.data.baranyi %>% filter(FALSE)
baranyi.curves.2 = baranyi.curves %>% filter(FALSE)
for (this.real.lag in real.lags.baranyi) {
  baranyi.simulated.data.basic = data.frame(time = times,
                                            biomass = Baranyi.Solution(times, LOG10Nmax, growth.rate.baranyi, LOG10N0, this.real.lag))
  
  lag.df.obj = suppressMessages(Get.Lag.Fitting.Data.For.Noisy.Simulations(baranyi.simulated.data.basic,
                                                                           sd_range = sd_range,
                                                                           biomass.increase.threshold,
                                                                           Num.obs = Num.obs,
                                                                           MIN.BIOMASS= MIN.BIOMASS))
  
  baranyi.curves.2 = rbind(baranyi.curves.2,
                            lag.df.obj$curves %>% mutate(growth.rate = growth.rate.baranyi,
                                                         real.lag = this.real.lag,
                                                         time.interval = time.interval))
  all.lag.data.baranyi.2 = rbind(all.lag.data.baranyi.2,
                         lag.df.obj$lag.df %>% mutate(growth.rate = growth.rate.baranyi,
                                                      real.lag = this.real.lag,
                                                      time.interval = time.interval))
}


all.lag.data.baranyi.3 = all.lag.data.baranyi %>% filter(FALSE)
baranyi.curves.3 = baranyi.curves %>% filter(FALSE)
for (this.growth.rate in growth.rates.baranyi) {

  baranyi.simulated.data.basic = data.frame(time = times,
                                            biomass = Baranyi.Solution(times, LOG10Nmax, this.growth.rate, LOG10N0, real.lag.baranyi)) 
  
  lag.df.obj = suppressMessages(Get.Lag.Fitting.Data.For.Noisy.Simulations(baranyi.simulated.data.basic,
                                                                           sd_range = sd_range,
                                                                           biomass.increase.threshold,
                                                                           Num.obs = Num.obs,
                                                                           MIN.BIOMASS= MIN.BIOMASS))
  baranyi.curves.3 = rbind(baranyi.curves.3,
                            lag.df.obj$curves %>% mutate(growth.rate = this.growth.rate,
                                                         real.lag = real.lag.baranyi,
                                                         time.interval = time.interval))
  all.lag.data.baranyi.3 = rbind(all.lag.data.baranyi.3,
                         lag.df.obj$lag.df %>% mutate(growth.rate = this.growth.rate,
                                                      real.lag = real.lag.baranyi,
                                                      time.interval = time.interval))
}


all.lag.data.baranyi.4 = all.lag.data.baranyi %>% filter(FALSE)
baranyi.curves.4 = baranyi.curves %>% filter(FALSE)
for (this.time.interval in time.intervals) {
  baranyi.simulated.data.basic = data.frame(time = seq(0,Max.Time,this.time.interval),
                                            biomass = Baranyi.Solution(seq(0,Max.Time,this.time.interval), 
                                                                       LOG10Nmax, growth.rate.baranyi, LOG10N0, real.lag.baranyi))
  
  
  lag.df.obj = suppressMessages(Get.Lag.Fitting.Data.For.Noisy.Simulations(baranyi.simulated.data.basic,
                                                                           sd_range = sd_range,
                                                                           biomass.increase.threshold,
                                                                           Num.obs = Num.obs,
                                                                           MIN.BIOMASS= MIN.BIOMASS))
  baranyi.curves.4 = rbind(baranyi.curves.4,
                            lag.df.obj$curves %>% mutate(growth.rate = growth.rate.baranyi,
                                                         real.lag = real.lag.baranyi,
                                                         time.interval = this.time.interval))
  all.lag.data.baranyi.4 = rbind(all.lag.data.baranyi.4, lag.df.obj$lag.df %>% mutate(growth.rate = growth.rate.baranyi,
                                                                      real.lag = real.lag.baranyi,
                                                                      time.interval = this.time.interval))
}

all.lag.data.baranyi.2 = all.lag.data.baranyi.2 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.baranyi.3 = all.lag.data.baranyi.3 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.baranyi.4 = all.lag.data.baranyi.4 %>% mutate(obs.minus.real.lag = lag - real.lag)

saveRDS(all.lag.data.baranyi.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.baranyi.varying.lags.rds"))
saveRDS(all.lag.data.baranyi.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.baranyi.varying.growth.rate.rds"))
saveRDS(all.lag.data.baranyi.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.baranyi.varying.time.interval.rds"))
saveRDS(baranyi.curves.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.baranyi.varying.lags.curves.rds"))
saveRDS(baranyi.curves.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.baranyi.varying.growth.rate.curves.rds"))
saveRDS(baranyi.curves.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.baranyi.varying.time.interval.curves.rds"))




