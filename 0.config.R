rm(list = ls())
set.seed(1)
library(data.table)
library(ggplot2)
library(deSolve)
library(DEoptim)
library(dplyr)
library(miLAG)
library(tidyr)
library(scales)


# Monika:
#THIS_PROJECT_PATH = "C:/Users/monik/Documents/_PROJEKTY i LAB/_2024 - LAG - PART 2/"
#OUTPUT_FIGS_PATH =sprintf("%soutput/2024.08.28 - figs and data/", THIS_PROJECT_PATH) 
#source(sprintf("%scode/in_silico analysis/lags_helper.R", THIS_PROJECT_PATH))
# Bogna
THIS_PROJECT_PATH = "/Users/bsmug/Library/CloudStorage/OneDrive-UniwersytetJagielloński/Lags_part2/"

OUTPUT_FIGS_PATH =sprintf("%soutput/Figures_2025_06_24/", THIS_PROJECT_PATH) 
EXPERIMENTAL.DATA.PATH = sprintf("%sexperimental data/", THIS_PROJECT_PATH)
NOLAG.DATA.PATH = sprintf("%sexperimental data/BY_a_lag_data.txt", THIS_PROJECT_PATH)

source(sprintf("%scode/lag_calc_methods_comparison/lags_helper.R", THIS_PROJECT_PATH))
dir.create(OUTPUT_FIGS_PATH)

# GLOBAL PARAMS
biomass.increase.threshold = 10^3
Max.Time = 24
method.labels = c(
  "biomass increase",
  "max growth acceleration",
  "tangent to max growth point",
  "tangent to max growth line",
  "par. fitting to logistic model",
  "par. fitting to baranyi model"
)