# Illustrative example based on Wyman et al., JCCP (2008)

# library(tidyverse)
library(lme4)
library(boot)
library(bootmlm)
source("compute_mlm_d.R")

dat <- read.csv("illustrative_revise.dat")

m1 <- lmer(y ~ treat + (1 | gp_id), data = dat)
summary(m1)

set.seed(171040)  # for replicable results
# ANOVA-based
ComputeES(dat, "y", "treat", "gp_id", type = "ss", ci = "central")
# Model-based
ComputeES(dat, "y", "treat", "gp_id", type = "model", ci = "central")
# Parametric bootstrap
ComputeES(dat, "y", "treat", "gp_id", type = "model",
          boot = "parametric", ci = "all",
          nsim = 1999L)
# Residual bootstrap
ComputeES(dat, "y", "treat", "gp_id", type = "model",
          boot = "residual", ci = "all",
          nsim = 1999L)
# Residual (CGR) bootstrap
ComputeES(dat, "y", "treat", "gp_id", type = "model",
          boot = "residual_cgr", ci = "all",
          nsim = 1999L)
# Case bootstrap
ComputeES(dat, "y", "treat", "gp_id", type = "model",
          boot = "case", ci = "all",
          nsim = 1999L)
# Ignore multilevel structure
psych::cohen.d(dat[c("y", "treat")], "treat")
