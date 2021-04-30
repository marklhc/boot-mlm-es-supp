ComputeICC <- function(df, y_var, treat_var, id_var) {
  # Function to compute intraclass correlation using the ANOVA approach.
  # df = input data frame, 
  # y_var = string for the name of the outcome variable
  # treat_var = string for the name of the treatment variable
  # id_var = string for the name of the group membership variable
  N <- nrow(df)
  J <- length(unique(df[[id_var]]))
  form <- paste0(y_var, " ~ ", treat_var, " + Error(", id_var, ")")
  sumfit <- summary(aov(as.formula(form), data = df))
  msb <- sumfit[[1]][[1]]["Residuals", "Mean Sq"]
  msw <- sumfit[[2]][[1]]["Residuals", "Mean Sq"]
  n <- N / J
  icc <- (msb - msw) / (msb + (n - 1) * msw)
  if (icc < 0) icc <- 0
  icc
}

vcov_sig <- function(x) {
  # Function to obtain asymptoic variance of sigma^2, the level-1 term
  # x = fitted object by lme4::lmer()
  dd <- bootmlm::devfun_mer(x)
  n_th <- length(x@theta)
  qs <- lengths(x@cnms)
  dd2 <- function(vc) {
    sigma <- sqrt(vc[n_th + 1])
    th_sig <- c(lme4::Vv_to_Cv(vc, n = qs, s = sigma), sigma)
    dd(th_sig)
  }
  dd3 <- function(s2) {
    dd2(c(0, s2))
  }
  hess <- numDeriv::hessian(dd3, sigma(x)^2)
  vv <- 2 / hess
  vv
}

mle_es <- function(model, treat_var, vd = FALSE) {
  # Model-based effect size estimate and sampling variance estimate 
  # (if vd = TRUE)
  sigma2T <- (model@theta^2 + 1) * sigma(model)^2
  gam10 <- fixef(model)[[treat_var]]
  est_dT <- gam10 / sqrt(sigma2T)
  if (vd) {
    v_gam10 <- vcov(model)[treat_var, treat_var]
    v_s2T <- try(
      sum(bootmlm::vcov_vc(model, sd_cor = FALSE, print_names = FALSE))
    )
    if (inherits(v_s2T, "try-error")) {
      v_s2T <- vcov_sig(model)
    }
    const_c <- v_s2T / 2 / sigma2T^2
    var_dT <- v_gam10 / sigma2T + est_dT^2 * const_c / 2
    return(c(est_dT, var_dT))
  }
  est_dT
}

extract_bootci <- function(boo_ci) {
  unlist(lapply(boo_ci[-c(1:3)], function(mat) {
    out <- tail(mat[1, ], 2L)
    names(out) <- c("ll", "ul")
    out
  }))
}

ComputeES <- function(df, y_var, treat_var, id_var, model = NULL,
                      type = c("ss", "model"), 
                      boot = "none", 
                      nsim = 1L,
                      ci = "central", 
                      L = NULL, ...) {
  # Function to compute multilevel effect size.
  # df = input data frame, 
  # y_var = string for the name of the outcome variable
  # treat_var = string for the name of the treatment variable
  # id_var = string for the name of the group membership variable
  # model = fitted object from lme4::lmer(). If NULL, a random-intercept model
  #         will be run first on the data with treat_var as the predictor
  # type = ANOVA-based or model-based CI. Use "model" for all bootstrap methods
  # boot = type of bootstrap to be passed on to bootmlm::bootstrap_mer()
  # nsim = number of bootstrapped samples
  # ci = type of confidence intervals, to be passed on to
  #      bootmlm::bootstrap_mer()
  # L = empirical influence function to be used for BCa CI. If NULL, L will be
  #     estimated using the grouped jackknife. 
  type <- match.arg(type)
  
  if (type == "ss") {
    N_T <- sum(df[[treat_var]] == 1)
    N_C <- sum(df[[treat_var]] == 0)
    N <- N_T + N_C
    J <- length(unique(df[[id_var]]))
    cs <- N / J
    form <- paste0(y_var, " ~ ", treat_var, " + Error(", id_var, ")")
    sumfit <- summary(aov(as.formula(form), data = df))
    msb <- sumfit[[2]][[1]][treat_var, "Mean Sq"]
    msw <- sumfit[[2]][[1]]["Residuals", "Mean Sq"]
    icc <- (msb - msw) / (msb + (cs - 1) * msw)
    if (icc < 0) icc <- 0
    
    const_b <- 1 - (2 * (cs - 1) * icc) / (N - 2)
    s2_T <- with(df, as.vector(crossprod(y - ave(y, treat))) / (N - 2))
    mean_diff <- with(df, as.vector(diff(tapply(y, treat, mean))))
    est_dT <- mean_diff * sqrt(const_b / s2_T)
    if ("central" %in% ci) {
      const_a <- 1 + (cs - 1) * icc
      const_c <- ((N - 2) * (1 - icc)^2 + cs * (N - 2 * cs) * icc^2 + 
                    2 * (N - 2 * cs) * icc * (1 - icc)) / (N - 2)^2
      N_tilde <- N_T * N_C / (N_T + N_C)
      var_dT <- const_a / N_tilde + const_c * est_dT^2 / 2 / const_b^2
    }
  } else if (type == "model") {
    if (missing(model)) {
      form <- paste0(y_var, " ~ ", treat_var, " + (1 | ", id_var, ")")
      model <- lmer(as.formula(form), data = df)
    }
    if ("central" %in% ci) {
      dT_vdT <- mle_es(model, treat_var, vd = TRUE)
      est_dT <- dT_vdT[1]
      var_dT <- dT_vdT[2]
    } else {
      est_dT <- mle_es(model, treat_var, vd = FALSE)
    }
  }
  out <- c(d = est_dT)
  if ("central" %in% ci) {
    ci_dT <- est_dT + c(-1, 1) * qnorm(.975) * sqrt(var_dT)
    names(ci_dT) <- c("central_ll", "central_ul")
    out <- c(out, ci_dT)
  } 
  if (boot != "none" & type == "model") {
    ci_type <- setdiff(ci, "central")
    if ("stud" %in% ci_type | ci_type == "all") {
      d_fun <- function(x) mle_es(x, treat_var = treat_var, vd = TRUE)
    } else {
      d_fun <- function(x) mle_es(x, treat_var = treat_var, vd = FALSE)
    }
    boo_d <- bootmlm::bootstrap_mer(model, d_fun, nsim, type = boot, ...)
    if (missing(L) & ("bca" %in% ci_type | ci_type == "all")) {
      d_fun2 <- function(x) mle_es(x, treat_var = treat_var, vd = FALSE)
      L <- bootmlm::empinf_mer(model, d_fun2, index = 1L)
    }
    boo_ci <- boot.ci(boo_d, type = ci_type, L = L)
    ci_dT <- extract_bootci(boo_ci)
    d_boot <- 2 * est_dT - mean(boo_d$t[ , 1], na.rm = TRUE)
    out <- c(out, d_boot = d_boot, ci_dT)
  }
  out
}
