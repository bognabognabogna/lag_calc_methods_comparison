Get.Theme = function(text.size = 12) {
  my_theme = theme(
    panel.grid.major = element_blank(),#element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_blank(),#element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text.size),
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black")
  )
  return(my_theme)
}



#' get_init_pars_baranyi: copied from milag
#'
#' Finds reasonable approximation for baranyi growth curve parameters (init_mumax, lag) based on the growth curve and some initial values
#' These approximations will be used as the initial values for the proper optimization algorithm run later.
#' @param data_this_curve data from one specific growth curve with these two columns: time and biomass
#' @param this_n0 the initial biomass
#' @param init_lag initial value for the lag parameter
#' @param init_gr_rate initial value for the growth rate
#' @param min_b defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(mina TO minb)
#' @param min_a defaults to 0.8
#' @returns list of parameters: init_mumax, init_lag
get_init_pars_baranyi <- function(data_this_curve, this_n0, init_lag, init_gr_rate, min_b = 0.2, min_a = 0.8) {
  if (is.null(init_lag)) {
    init_lag <- calc_lag(data_this_curve, method = "tangent", pars = get_def_pars()) %>% pull(lag) %>% unique() %>% as.numeric()
  }
  
  if (is.null(init_gr_rate)) {
    data_this_curve_exp <- data_this_curve %>%
      mutate(
        max_biomass = max(biomass),
        min_threshold = this_n0 + min_b * (max_biomass - this_n0),
        max_threshold = this_n0 + min_a * (max_biomass - this_n0)) %>%
      # take only the points that are in between min and max
      filter(biomass <= max_threshold & biomass >= min_threshold)
    data_this_curve_exp$logdata <- log(data_this_curve_exp$biomass/this_n0)
    if (nrow(data_this_curve_exp %>% filter(!is.na(time) & !is.na(logdata))  > 0)) {
      mod <- lm(logdata ~ time, data = data_this_curve_exp)
      # this growth rate is assuming an exponential model so it will be generally underestimated
      init_gr_rate <- mod$coefficients[2] %>% unname()
      # we have real r = r(1-N/K) so let us take
      n_mid <- median(data_this_curve_exp$biomass)
      init_mumax <- init_gr_rate
    } else {
      init_mumax <- 0.1
    }
  } else {
    init_mumax <- init_gr_rate
  }
  return(list(init_mumax = init_mumax, init_lag = init_lag))
}


Baranyi.Solution = function(t, LOG10Nmax, mumax, LOG10N0, lag) {
  LOG10N = LOG10Nmax + log10((-1 + exp(mumax*lag) + exp(mumax*t))/(exp(mumax*t) - 1 + exp(mumax*lag) * 10^(LOG10Nmax - LOG10N0)))
  N = 10^LOG10N
  return(N)
}

Fit.To.Baranyi = function(data_this_curve, max_iter = 500, init_lag, init_gr_rate) {
  data_this_curve_for_model <- data_this_curve %>%
    mutate(LOG10N = log10(biomass), t = time) %>%
    select(LOG10N, t)
  n0 = data_this_curve %>% arrange(time) %>% head(1)
  init_LOG10N0 <- log10(n0$biomass)
  init_LOG10Nmax <- max(data_this_curve_for_model$LOG10N)
  init_pars_baranyi <- get_init_pars_baranyi(data_this_curve, n0$biomass, init_lag, init_gr_rate)
  nlsres_LM <- nlsLM(formula = baranyi,
                                   data = data_this_curve_for_model,
                                   start = list(lag=init_pars_baranyi$init_lag, 
                                                mumax= init_pars_baranyi$init_mumax, 
                                                LOG10N0 = init_LOG10N0, 
                                                LOG10Nmax = init_LOG10Nmax),
                                   control = nls.control(maxiter = max_iter),
                                   lower = c(0,0,0, 0))
  return(nlsres_LM)
}


sumLeastSquaresFitLogisticToDeoptim = function(param, N0, data) {
  growth.rate = param[1]
  K = param[2]
  lag = param[3]
  simulatedN = Simulate.Logistic.With.Lag(N0, growth.rate, K, lag, data$time)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}

sumLeastSquaresFitMonodToDeoptim = function(param, N0, G0, data) {
  a = param[1]
  Vh = param[2]
  Kh = param[3]
  lag = param[4]
  simulatedN = Simulate.Monod.With.Lag(a,Vh,Kh, N0, G0, lag, data$time) 
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}


Fit.To.Logistic.With.Lag = function(N0, data, lower.bound = c(0.01,1,0), upper.bound = c(2,50,10)) {
  deopticontrol = DEoptim.control(itermax = 100, reltol = 10^(-8), trace = 100)
  deoptim_out = DEoptim(fn = function(param) {sumLeastSquaresFitLogisticToDeoptim(param, N0, data)},
                        lower = lower.bound,
                        upper = upper.bound,
                        control = deopticontrol)
  return(deoptim_out)
}

Fit.To.Monod.With.Lag = function(N0, G0, data, lower.bound = c(0.01,1,1, 0), upper.bound = c(2,50,50, 10),
                                 reltol = 10^(-8), itermax = 100) {
  deopticontrol = DEoptim.control(itermax = itermax, reltol = reltol, trace = 100)
  deoptim_out = DEoptim(fn = function(param) {sumLeastSquaresFitMonodToDeoptim(param, N0, G0, data)},
                        lower = lower.bound,
                        upper = upper.bound,
                        control = deopticontrol)
  return(deoptim_out)
}


###################### Data Simulation ##############################

logistic_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ndot   = r*N*(1-N/K)
    return(list(c(Ndot)))
  })
}

Simulate.Logistic.With.Lag = function(N0, growth.rate, K, lag, times) {
  #logistic model solution derived for example here
  #https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/08%3A_Introduction_to_Differential_Equations/8.4%3A_The_Logistic_Equation
  N = (times >= lag)*(N0*K*exp(growth.rate*(times-lag))/(K - N0 + N0*exp(growth.rate*(times - lag)))) +
      (times < lag)*N0
  logistic.simulated.data = data.frame(time = times, biomass = N)
  return(logistic.simulated.data)
}

Simulate.Logistic.With.Lag.Old = function(N0, growth.rate, K, lag, times) {
  pars <- c(r = growth.rate, K=K)
  inits = c(N = N0)
  lag.data = data.frame(time = times) %>% mutate(biomass = N0) %>% filter(time <= lag)
  lagged.times = times - lag
  positive.lagged.times = lagged.times[lagged.times >= 0]
  logistic.simulated.data <- as.data.frame(ode(inits, positive.lagged.times, logistic_model, pars)) %>%
    select(time = time, biomass = N) %>%
    mutate(time = time + lag) %>%
    filter(time > lag & time <= max(times))
  logistic.simulated.data = rbind(lag.data, logistic.simulated.data)
  return(logistic.simulated.data)
}


baranyi_and_roberts_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    f = 1-(N/K)
    alpha = Q/(1+Q)
    Qdot = v*Q
    Ndot   = mu_max*alpha*f*N
    return(list(c(Ndot, Qdot)))
  })
}


Monod_bacteria_growth_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G)

    if (Time <= lag) {
      Ndot = 0
      Gdot = 0
    } else {
      Ndot = a*Jg*N
      Gdot   = -Jg*N
    }

    return(list(c(Gdot, Ndot)))
  })
}

Simulate.Monod.With.Lag = function(a,Vh,Kh, N0, G0, lag,  times) {
  pars <- c(a = a, Vh = Vh, Kh = Kh, lag = lag)
  inits = c(G=G0, N=N0)
  lag.data = data.frame(time = times) %>% mutate(biomass = N0) #%>% filter(time <= lag)

  monod.simulated.data <- as.data.frame(ode(inits, times, Monod_bacteria_growth_model, pars)) %>%
    select(time = time, biomass = N) #%>%
    #mutate(time = time + lag) %>%
    #filter(time > lag)
  #monod.simulated.data = rbind(lag.data, monod.simulated.data)
  return(monod.simulated.data)
}


Simulate.Exponential.Data.With.Lag = function(N0, growth.rate, lag, times) {
  exponential.growth.simulated.data = data.frame(
    time = times) %>%
    mutate(lagged.time = if_else(time - lag < 0, 0, time - lag)) %>%
    mutate(biomass = N0*exp(growth.rate*lagged.time))
  return(exponential.growth.simulated.data)
}



Plot.Lag.Fit = function(data.new, print.lag.info = TRUE) {
  data.new = data.new %>%
    group_by(curve_id) %>%
    mutate(x.mid = mean(time),
           lag.info = paste0("Lag = ", round(lag, 3), " [h]."),
           log.biomass = log(biomass),
           log.10.tangent.point = log10(tangent_point),
           log.10.biomass = log10(biomass),
           log.10.predicted = log10(predicted_data),
           log.10.threshold = log10(threshold),
           y.max.for.curve = max(log.10.biomass),
           y.min.for.curve = min(log.10.biomass),
           log10N0 = log10(exp(log(n0))),
           text.y = 1.005*y.max.for.curve,
           max.second.deriv.b = max(second_deriv_b, na.rm = TRUE),
           min.second.deriv.b = min(second_deriv_b, na.rm = TRUE),
           second.deriv.b.scaled = (second_deriv_b - min.second.deriv.b)/(max.second.deriv.b - min.second.deriv.b)*(y.max.for.curve - y.min.for.curve) + y.min.for.curve
           #y.limit = 1.1*y.max.for.curve
           ) %>%
    ungroup() %>%
    mutate(min.log10N0 = min(log.10.biomass),
           max.log10N0 = max(log.10.biomass),
           log10.intercept = line_intercept/log(10),
           log10.slope = line_slope/log(10))

  max.time = max(data.new$time)
    #mutate(curve_id = paste0(curve_id, ":\n", lag.info))

  #coef.diff = max(data.new$y.max, na.rm = TRUE)/max(data.new$diff, na.rm = TRUE)
  #coef.second.deriv.b = max(data.new$max.log10N0, na.rm = TRUE)/max(data.new$second.deriv.b, na.rm = TRUE)

  size.N0.line = 1
  size.lag.line = 1
  g = ggplot(data.new)  +
    geom_vline(aes(xintercept = lag), size = size.lag.line, col = "red", linetype = "dashed") +
    geom_line(aes(x= time, y = log.10.biomass), col = "blue") +
    #geom_point(aes(x= time, y = log10.biomass), col = "blue") +
    geom_point(aes(x= time, y = log.10.tangent.point), col = "darkgreen", size = 2) +
    geom_line(aes(x= time, y = log.10.predicted), col = "darkgreen") +
    geom_line(aes(x= time, y = log.10.threshold), col = "darkgreen") +
    geom_line(aes(x=time, y = second.deriv.b.scaled), col = "darkgreen", alpha = 0.5) +
    geom_hline(aes(yintercept = log10N0), size = size.N0.line, col = "black") +
    geom_abline(aes(intercept = log10.intercept, slope = log10.slope), col = "darkgreen") +
    xlab("time [h]") +
    xlim(c(0, max.time)) +
    ylab("Log10(biomass)") +
    #ylim(c(min(data.new$log.10.biomass), max(data.new$log.10.biomass))) +
    facet_grid(curve_id~lag_calculation_method, scales = "free") +
  theme(axis.text.y.right=element_text(colour="black"),
          axis.text.y=element_text(colour="blue"),
          axis.title.y=element_text(colour="blue"),
          axis.title.y.right=element_text(colour="black"))

  if (print.lag.info) {
    g = g +  geom_text(aes(x=x.mid, y = text.y, label = lag.info), size = 6, col = "red")
  }
  #g =  g +
  #  geom_line(aes(x = time, y = diff), alpha = 0.5, col = "black")  +
  #  geom_line(aes(x = time, y = second.deriv.b), alpha = 0.5, col = "black")

 # if ("second.deriv.b" %in% colnames(data.new)) {
  #  g= g +
  #  geom_line(aes(x = time, y = second.deriv.b)) +
  #  scale_y_continuous(name = "Log(biomass)",
  #                     sec.axis = sec_axis(~.,
  #                                         name = "Second derivative of log(biomass)"))
  #
  #}

  #else if ("diff" %in% colnames(data.new)) {
  #  g =  g + geom_line(aes(x = time, y = diff))  +
  #    scale_y_continuous(name = "Log(biomass)",
  #                       sec.axis = sec_axis(~.,
  #                                          name = "Increase in log(biomass)"))
  # }

  return(g)
}





Get.Lag.Fitting.Data.For.Noisy.Simulations = function(simulated.data,
                                                      sd_range = seq(0.0, 0.5, 0.05),
                                                      biomass.increase.threshold,
                                                      Num.obs = 100,
                                                      curve_name = "curve",
                                                      MIN.BIOMASS) {
  curves = data.frame(time = numeric(0), biomass = numeric(0), curve_id = character(0), sd = numeric(0))
  lag.df = data.frame(curve_id = character(0),
                      lag = numeric(0),
                      lag_calculation_method = character(0),
                      sd = numeric(0))
  N0 = simulated.data$biomass[1]
  for (this.sd in sd_range) {
    curves.this.sd = data.frame(time = numeric(0), biomass = numeric(0), curve_id = character(0))
    for (i in 1:Num.obs) {
      noise = rnorm(n = nrow(simulated.data), mean = 0, sd = this.sd*N0)
      curve_i = simulated.data %>%
        mutate(biomass = biomass + noise,
               curve_id = i,
               sd = this.sd) %>%
        as.data.table()
      # MAKE SURE THE SIMULATED DATA IS >= 0
      # if biomass is exactly zero we have issues with LOG and start getting Inf values. Biomass needs to be bigger than 0
      # DON"T use na.rm = T because that will convert NA to 0!
       curve_i[, biomass := ifelse(is.na(biomass), biomass, pmax(biomass, MIN.BIOMASS))]
       curve_i = curve_i %>% as.data.frame()
      curves.this.sd = rbind(curves.this.sd, curve_i)
    }
    data.all.with.lag = get_all_methods_lag(curves.this.sd, biomass.increase.threshold)
    lag.data = data.all.with.lag %>% distinct(curve_id, lag, lag_calculation_method)
    lag.df = rbind(lag.df, lag.data %>% mutate(sd = this.sd))
    curves = rbind(curves, curves.this.sd %>% mutate(sd = this.sd))
  }
  return(list(lag.df = lag.df, curves = curves))
}


