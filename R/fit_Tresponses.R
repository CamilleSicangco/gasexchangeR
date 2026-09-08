# Temperature-response formula builders ---------------------------------------
#
# Adapted from scripts by Belinda Medlyn, Dushan Kumarathunge, and John Drake.
# These two helpers replace 11 near-duplicate top-level formula objects
# (fVc, fVJ, fVc.a, fVJ.a, f_topt_V, f_topt_J, fA, fA250, fA275, fAc250,
# fAj250) that differed only in which response-variable name appeared on the
# left-hand side.

#' Build an Arrhenius-family temperature-response formula
#'
#' @param response Name of the response-variable column, e.g. \code{"Vcmax"}
#'   or \code{"Jmax"}
#' @param form One of \code{"standard"} (plain Arrhenius), \code{"delS_peaked"}
#'   (peaked Arrhenius with an entropy term), or \code{"topt_peaked"} (peaked
#'   form parameterized directly by \code{Topt})
#'
#' @return A \code{formula}
#' @export
arrhenius_formula = function(
  response,
  form = c("standard", "delS_peaked", "topt_peaked")
)
{
  form = match.arg(form)

  rhs = switch(form,
    standard = "k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK))",
    delS_peaked = paste0(
      "k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * ",
      "(1+exp((298.15*delS - 200)/(298.15*0.008314))) / ",
      "(1+exp((TsK*delS-200)/(TsK*0.008314)))"
    ),
    topt_peaked = paste0(
      "Ropt*((200*exp(Ea*(TsK-Topt)/(TsK*0.008314*Topt))/",
      "(200-(Ea*(1-exp(200*(TsK-Topt)/(TsK*0.008314*Topt)))))))"
    )
  )

  as.formula(paste(response, "~", rhs))
}

#' Build a quadratic-Topt temperature-response formula
#'
#' @param response Name of the response-variable column, e.g. \code{"A250"},
#'   \code{"Ac250"}, \code{"Aj250"}, \code{"A275"}, or \code{"A"}
#'
#' @return A \code{formula}
#' @export
quadratic_topt_formula = function(response)
{
  as.formula(paste0(response, " ~ Aopt-(b*(Tleaf-Topt)^2)"))
}

# nlme random-effect helpers ---------------------------------------------------
#
# fit_temperature_response(), fit_arrhenius(), and fitAvT() all fit a model
# with a per-group random intercept via nlme(), grouped by a column whose
# name is now caller-specified (random_effect_col) rather than hardcoded.
# nlme can't estimate a random effect from a single group, so these helpers
# let each function fall back to an equivalent fixed-effects-only nls() fit
# rather than erroring when random_effect_col only has one level.

#' Whether a column has more than one distinct value
#' @noRd
has_multiple_levels = function(data, col)
{
  length(unique(data[[col]])) >= 2
}

#' Build a "~ 1 | group" random-intercept formula for nlme
#' @noRd
random_intercept_formula = function(lhs, group_col)
{
  as.formula(paste0(lhs, " ~ 1 | ", group_col))
}

#' Warn that a random effect is being dropped for lack of grouping levels
#' @noRd
warn_no_random_effect = function(group_col)
{
  message(
    "Only one level of '", group_col, "' present in the data -- fitting ",
    "without a random effect on ", group_col, "."
  )
}

#' Fixed-effect coefficient table from either an nlme or nls fit
#'
#' \code{summary(nlme_fit)$tTable} and \code{summary(nls_fit)$coefficients}
#' are both parameter-by-statistic matrices with \code{Value}/\code{Estimate}
#' as column 1 and \code{Std.Error} as column 2 (nlme has extra columns
#' after that; nls doesn't) -- so positional indexing into the first two
#' columns behaves identically either way.
#' @noRd
get_coef_table = function(fit)
{
  if (inherits(fit, "nlme")) summary(fit)$tTable else summary(fit)$coefficients
}

#' p-value from anova() on a pair of nlme or nls fits
#'
#' nlme::anova.lme() and stats::anova.nls() name their p-value column
#' differently ("p-value" vs. "Pr(>F)"), so look it up by name rather than
#' position.
#' @noRd
get_anova_pvalue = function(an)
{
  pcol = if ("p-value" %in% names(an)) "p-value" else "Pr(>F)"
  an[[pcol]][2]
}

# Fit Vcmax/Jmax temperature responses ----------------------------------------

#' Fit an Arrhenius/peaked-Arrhenius temperature response
#'
#' Fits a peaked Arrhenius model (with an entropy term) and a standard
#' Arrhenius model to \code{Vcmax} or \code{Jmax} vs. leaf temperature, with
#' a per-\code{random_effect_col} random effect, and compares the two via
#' AIC/anova. If \code{random_effect_col} only has one level in \code{dat},
#' fits without a random effect (via \code{nls}) instead of erroring.
#' Replaces the previously separate \code{fitpeaked()} (Vcmax) and
#' \code{fitpeakedJ()} (Jmax) functions.
#'
#' @param dat A data frame with columns \code{TsK} (leaf temperature, K),
#'   \code{random_effect_col}, and \code{response}
#' @param response Which response variable to fit: \code{"Vcmax"} or
#'   \code{"Jmax"}
#' @param return Which model's parameters to return: \code{"Arr"} (standard
#'   Arrhenius) or \code{"Peak"} (peaked Arrhenius, plus its fitted Topt)
#' @param start Starting values for the peaked-Arrhenius fit
#' @param random_effect_col Name of the column identifying groups for the
#'   per-group random effect (e.g. individual/curve identity)
#'
#' @return A one-row data frame of fitted parameters and standard errors,
#'   plus \code{R2_Arr}, \code{R2_Peak}, \code{AIC_Arr}, \code{AIC_Peak}, and
#'   the anova p-value comparing the two models
#' @export
fit_temperature_response = function(
  dat,
  response = c("Vcmax", "Jmax"),
  return = c("Arr", "Peak"),
  start = c(k25 = 100, Ea = 60, delS = 0.64),
  random_effect_col = "curve_ID"
)
{
  response = match.arg(response)
  return = match.arg(return)
  suffix = substr(response, 1, 1)

  dat = dat |> filter(!is.na(!!sym(response)))

  use_random_effect = has_multiple_levels(dat, random_effect_col)
  if (!use_random_effect) warn_no_random_effect(random_effect_col)

  # Peaked Arrhenius model
  Vc <- if (use_random_effect) {
    nlme(model = arrhenius_formula(response, "delS_peaked"),
         fixed = (k25 + Ea + delS ~ 1),
         random = random_intercept_formula("k25", random_effect_col),
         data = dat,
         start = start)
  } else {
    nls(arrhenius_formula(response, "delS_peaked"), data = dat, start = start)
  }

  res1 <- data.frame(t(get_coef_table(Vc)[1:6]))
  names(res1) <- c(
    paste0(response, "25"), paste0("Ea", suffix), paste0("dels", suffix),
    paste0(response, "25.se"), paste0("Ea", suffix, ".se"), paste0("dels", suffix, ".se")
  )

  # Standard Arrhenius model
  Vr <- if (use_random_effect) {
    nlme(model = arrhenius_formula(response, "standard"),
         fixed = (k25 + Ea ~ 1),
         random = random_intercept_formula("k25", random_effect_col),
         data = dat,
         start = c(k25 = 100, Ea = 40))
  } else {
    nls(arrhenius_formula(response, "standard"), data = dat, start = c(k25 = 100, Ea = 40))
  }

  res <- data.frame(t(get_coef_table(Vr)[1:4]))
  names(res) <- c(
    paste0(response, "25"), paste0("Ea", suffix),
    paste0(response, "25.se"), paste0("Ea", suffix, ".se")
  )

  an <- anova(Vr, Vc)
  AIC_Arr <- AIC(Vr)
  AIC_Peak <- AIC(Vc)
  prob <- get_anova_pvalue(an)

  r1 <- cor(fitted(Vc), dat[[response]])
  R2_Peak <- r1 * r1

  r2 <- cor(fitted(Vr), dat[[response]])
  R2_Arr <- r2 * r2

  # Test for normality of residuals (diagnostic only, not returned)
  shapiro.test(residuals(Vc))

  # Topt of the response, with SE (peaked form parameterized by Topt directly)
  VcTopt <- if (use_random_effect) {
    nlme(model = arrhenius_formula(response, "topt_peaked"),
         fixed = (Ropt + Ea + Topt ~ 1),
         random = random_intercept_formula("Ropt", random_effect_col),
         data = dat,
         start = c(Ropt = 150, Ea = 60, Topt = 310))
  } else {
    nls(arrhenius_formula(response, "topt_peaked"), data = dat,
        start = c(Ropt = 150, Ea = 60, Topt = 310))
  }

  VcTopt_coefs <- get_coef_table(VcTopt)
  res2 <- data.frame(cbind(VcTopt_coefs[1], VcTopt_coefs[3] - 273.15,
                           VcTopt_coefs[4], VcTopt_coefs[6]))
  names(res2) <- c(paste0("opt", suffix), paste0("Topt", suffix),
                   paste0("opt", suffix, ".se"), paste0("Topt", suffix, ".se"))

  if (return == "Peak") {
    param_peak <- cbind(res1, res2)
    param_peak$R2_Arr <- R2_Arr
    param_peak$R2_Peak <- R2_Peak
    param_peak$AIC_Arr <- AIC_Arr
    param_peak$AIC_Peak <- AIC_Peak
    param_peak[[paste0("prob_", suffix)]] <- prob
    return(param_peak)
  } else {
    param_arr <- res
    param_arr$R2_Arr <- R2_Arr
    param_arr$R2_Peak <- R2_Peak
    param_arr$AIC_Arr <- AIC_Arr
    param_arr$AIC_Peak <- AIC_Peak
    param_arr[[paste0("prob_", suffix)]] <- prob
    return(param_arr)
  }
}

#' Fit a standard Arrhenius temperature response only
#'
#' Fits only the standard (non-peaked) Arrhenius model to \code{Vcmax} or
#' \code{Jmax} vs. leaf temperature, with a per-\code{random_effect_col}
#' random effect — simpler than
#' \code{\link{fit_temperature_response}}\code{(return = "Arr")}, which also
#' fits the peaked model for comparison. If \code{random_effect_col} only
#' has one level in \code{dat}, fits without a random effect (via
#' \code{nls}) instead of erroring. Replaces the previously separate
#' \code{fitarh()} (Vcmax) and \code{fitarhJ()} (Jmax) functions.
#'
#' @param dat A data frame with columns \code{TsK} (leaf temperature, K),
#'   \code{random_effect_col}, and \code{response}
#' @param response Which response variable to fit: \code{"Vcmax"} or
#'   \code{"Jmax"}
#' @param random_effect_col Name of the column identifying groups for the
#'   per-group random effect (e.g. individual/curve identity)
#'
#' @return A named numeric vector of fitted parameters, standard errors, and
#'   \code{R2_Arr}
#' @export
fit_arrhenius = function(
  dat,
  response = c("Vcmax", "Jmax"),
  random_effect_col = "curve_ID"
)
{
  response = match.arg(response)
  suffix = substr(response, 1, 1)

  dat = dat |> filter(!is.na(!!sym(response)))

  use_random_effect = has_multiple_levels(dat, random_effect_col)
  if (!use_random_effect) warn_no_random_effect(random_effect_col)

  Vr <- if (use_random_effect) {
    nlme(model = arrhenius_formula(response, "standard"),
         fixed = (k25 + Ea ~ 1),
         random = random_intercept_formula("k25", random_effect_col),
         data = dat,
         start = c(k25 = 100, Ea = 40))
  } else {
    nls(arrhenius_formula(response, "standard"), data = dat, start = c(k25 = 100, Ea = 40))
  }

  res <- get_coef_table(Vr)[1:4]
  names(res) <- c(
    paste0(response, "25"), paste0("Ea", suffix),
    paste0(response, "25.se"), paste0("Ea", suffix, ".se")
  )

  r <- cor(fitted(Vr), dat[[response]])
  R2 <- r * r

  # Test for normality of residuals (diagnostic only, not returned)
  shapiro.test(residuals(Vr))

  param <- c(res, R2_Arr = R2)
  return(param)
}

# Fit A vs. leaf temperature (Topt model) -------------------------------------

#' Fit a quadratic Topt model to assimilation vs. leaf temperature
#'
#' Fits a quadratic Topt model (\code{A ~ Aopt - b*(Tleaf - Topt)^2}) to net
#' assimilation vs. Tleaf for a chosen fixed-Ci response variable, with or
#' without a per-\code{random_effect_col} random effect. If
#' \code{incl_random_effect} is \code{TRUE} but \code{random_effect_col}
#' only has one level in \code{df}, fits without a random effect instead of
#' erroring.
#'
#' @param df A data frame with columns \code{Tleaf}, \code{random_effect_col},
#'   and \code{var}
#' @param var Which response variable to fit: \code{"A250"}, \code{"Ac250"},
#'   \code{"Aj250"}, \code{"A275"}, or the fully dynamic \code{"A"}
#' @param incl_random_effect Whether to fit a per-\code{random_effect_col}
#'   random effect on \code{Aopt} (via \code{nlme}) or a fixed-effects-only
#'   fit (via \code{nls})
#' @param random_effect_col Name of the column identifying groups for the
#'   per-group random effect (e.g. individual/curve identity)
#'
#' @return A named numeric vector of fitted parameters, standard errors, and
#'   the temperature span (\code{Tspan}) for a target fraction of Aopt (see
#'   \code{\link{solve_Tspan}})
#' @export
fitAvT = function(
  df,
  var = c("A250", "Ac250", "Aj250", "A275", "A"),
  incl_random_effect = TRUE,
  random_effect_col = "curve_ID"
)
{
  var = match.arg(var)
  df = df |> filter(!is.na(!!sym(var)))

  formula_v = quadratic_topt_formula(var)
  start_vals = c(Aopt = max(df[[var]]), Topt = 25, b = 0.05)

  use_random_effect = isTRUE(incl_random_effect) && has_multiple_levels(df, random_effect_col)
  if (isTRUE(incl_random_effect) && !use_random_effect) warn_no_random_effect(random_effect_col)

  if (!use_random_effect) {
    A_Topt <- nls(formula_v, data = df, start = as.list(start_vals))
  } else {
    A_Topt <- nlme(model = formula_v,
                   fixed = (Aopt + b + Topt ~ 1),
                   random = random_intercept_formula("Aopt", random_effect_col),
                   data = df,
                   start = start_vals)
  }

  A_Topt2 <- summary(A_Topt)

  if (!use_random_effect) {
    Tspan <- solve_Tspan(Aopt = A_Topt2$coefficients[1],
                         Topt = A_Topt2$coefficients[2],
                         b = A_Topt2$coefficients[3])
    Tspan.se <- Tspan * sqrt((A_Topt2$coefficients[4] / A_Topt2$coefficients[1])**2 +
                               (A_Topt2$coefficients[6] / A_Topt2$coefficients[3])**2)
    results <- c(A_Topt2$coefficients[1:6], Tspan, Tspan.se)
  } else {
    Tspan <- solve_Tspan(Aopt = A_Topt2$tTable[1],
                         Topt = A_Topt2$tTable[3],
                         b = A_Topt2$tTable[2])
    Tspan.se <- Tspan * sqrt((A_Topt2$tTable[4] / A_Topt2$tTable[1])**2 +
                               (A_Topt2$tTable[5] / A_Topt2$tTable[2])**2)
    results <- c(A_Topt2$tTable[c(1, 3, 2, 4, 6, 5)], Tspan, Tspan.se)
  }

  # var -> (name prefix, "b"/"Tspan" suffix) — the only thing that varied
  # across the branches of the original var-specific implementations
  prefix = switch(var, A250 = "A250", Ac250 = "Ac", Aj250 = "Aj", A275 = "A275", A = "A")
  suffix = switch(var, A250 = "_250", Ac250 = "_c", Aj250 = "_j", A275 = "_275", A = "")

  names(results)[1:8] <- c(
    paste0(prefix, "opt"), paste0("Topt", prefix), paste0("b", suffix),
    paste0(prefix, "opt.se"), paste0("Topt", prefix, ".se"), paste0("b", suffix, ".se"),
    paste0("Tspan", suffix), paste0("Tspan", suffix, ".se")
  )

  return(results)
}

#' Fit an Arrhenius model to respiration vs. leaf temperature
#'
#' Adapted from Drake et al. 2015.
#'
#' @param dataresp A data frame with columns \code{Rd} and \code{Tleaf} (deg C)
#'
#' @return A one-row data frame with \code{Rday_25}, \code{Rday_25_SE},
#'   \code{Ea_Rday}, \code{Ea_Rday_SE}, and \code{Q10}; all \code{NA} if the
#'   fit fails
#' @export
fitRvT_Ea <- function(dataresp)
{
  Rgas <- 8.314 # J per mol per K
  rday <- try(nls(Rd ~ R_Tref * exp((Ea * 1000 / (Rgas * 298.15)) * (1 - 298.15 / (Tleaf + 273.15))),
                  start = list(R_Tref = 2, Ea = 20), data = dataresp), silent = TRUE)

  Rday_25 <- try(summary(rday)$coefficients[[1]], silent = TRUE)
  Rday_25_SE <- try(summary(rday)$coefficients[[3]], silent = TRUE)
  Ea_Rday <- try(summary(rday)$coefficients[[2]], silent = TRUE)
  Ea_Rday_SE <- try(summary(rday)$coefficients[[4]], silent = TRUE)
  Q10 <- try(exp(10 * Ea_Rday * 1000 / (Rgas * (25 + 273.15)^2)), silent = TRUE)

  param <- try(data.frame(Rday_25, Rday_25_SE, Ea_Rday, Ea_Rday_SE, Q10), silent = TRUE)
  try(names(param) <- c("Rday_25", "Rday_25_SE", "Ea_Rday", "Ea_Rday_SE", "Q10"), silent = TRUE)

  na_s <- data.frame(NA, NA, NA, NA, NA)
  names(na_s) <- c("Rday_25", "Rday_25_SE", "Ea_Rday", "Ea_Rday_SE", "Q10")

  if (is.data.frame(param)) {
    return(param)
  } else {
    return(na_s)
  }
}

#' Fit a mixed-model quadratic-in-log Arrhenius model to respiration vs. leaf temperature
#'
#' Adapted from Heskel et al. 2016. Fits \code{log(var) ~ Tleaf + Tleaf^2}
#' with a per-\code{replicate} (here, \code{sample}) random effect on the
#' intercept via \code{\link[lme4]{lmer}}, then back-transforms to
#' \code{Rdark} at 25 deg C. Distinct from \code{\link{fitRvT_Ea}} (a
#' Drake et al. 2015 Arrhenius \code{nls} fit with no random effect) —
#' the two use different models and different required columns.
#'
#' @param df A data frame with columns \code{sample}, \code{Tleaf} (deg C),
#'   and \code{var}, with one row per respiration measurement
#' @param var Name of the respiration-rate column to fit: one of
#'   \code{"R"}, \code{"Rdark"}, \code{"Rday"}, or \code{"Rd"}
#'
#' @return A named numeric vector: \code{a}/\code{b}/\code{c} (fixed-effect
#'   coefficients and their standard errors), \code{Rdark25}/\code{Rdark25.se},
#'   \code{r2}, and \code{Tleaf_min}/\code{Tleaf_max}
#' @export
fitRvT <- function(df, var = c("R", "Rdark", "Rday", "Rd"))
{
  var = match.arg(var)

  df = df |>
    filter(Tleaf < 50) |>
    mutate(ln_R = log(!!sym(var)), Tleaf2 = Tleaf^2)

  try(
    mod <- lmer(
      ln_R ~ Tleaf + Tleaf2 + (1 | sample),
      df)
  )
  R.1 = summary(mod)
  results <- (R.1$coefficients[1:6])

  Rdark25 = exp(results[1] + results[2] * 25 + results[3] * 25^2)
  Rdark25.se = sqrt((results[4])^2 + (results[5] * 25)^2 + (results[3] * 25^2)^2) * Rdark25

  r <- cor(fitted(mod), df[[var]])
  r2 <- r*r

  Tleaf_min = min(df$Tleaf)
  Tleaf_max = max(df$Tleaf)

  results = c(results, Rdark25, Rdark25.se, r2, Tleaf_min, Tleaf_max)

  names(results)[1:11] <- c("a","b","c","a.se","b.se","c.se", "Rdark25", "Rdark25.se", "r2", "Tleaf_min", "Tleaf_max")

  return(results)
}

#' Temperature span for a target fraction of Aopt
#'
#' Solves the quadratic Topt model
#' (\code{A = Aopt - b*(Tleaf - Topt)^2}) for the width of the temperature
#' window over which \code{A} stays above \code{span * Aopt}.
#'
#' @param Aopt Fitted optimum assimilation rate
#' @param Topt Fitted optimum temperature
#' @param b Fitted curvature parameter
#' @param span Target fraction of \code{Aopt} defining the window (default
#'   0.8, i.e. the 80\% window)
#'
#' @return The temperature span, same units as \code{Topt}
#' @export
solve_Tspan = function(Aopt, Topt, b, span = 0.8)
{
  rt = sqrt((1 - span) * Aopt / b)
  T1 = Topt - rt
  T2 = Topt + rt
  Tspan = T2 - T1
  return(Tspan)
}

# Fitting multiple groups of curves at once -----------------------------------

#' Fit a temperature-response model to multiple groups of curves at once
#'
#' Splits \code{data} by \code{group_cols} and applies \code{fit_fun} to
#' each group in turn (forwarding \code{...}), combining the results into a
#' single data frame with the grouping columns attached. Works with any of
#' this package's single-group temperature-response fitters --
#' \code{\link{fitRvT}}, \code{\link{fitRvT_Ea}}, \code{\link{fit_arrhenius}},
#' \code{\link{fit_temperature_response}}, \code{\link{fitAvT}} -- or a
#' custom function of your own, since it never inspects \code{data}'s column
#' names itself; it only needs \code{group_cols} to exist. Written to
#' replace ad hoc copies of this same split/fit/tryCatch/recombine pattern
#' that were starting to get written out by hand in downstream projects, one
#' project at a time.
#'
#' A "group" here means a set of curves to be pooled into one fit (e.g. every
#' curve for one species at one elevation) -- it is not a single curve's own
#' identifier (e.g. \code{curve_ID}/\code{sample}/\code{replicate}), which
#' most of these fitting functions already expect as a column within each
#' group and use internally for their own per-curve random effect. Pass that
#' column name as one of \code{fit_fun}'s own arguments (the same way you
#' would calling \code{fit_fun} directly), not as part of \code{group_cols}.
#'
#' @param data A data frame containing every group's curves, plus the
#'   column(s) named in \code{group_cols}
#' @param group_cols Character vector of column name(s) defining a pooling
#'   group (e.g. \code{c("sp_code", "elevation")}) -- everything sharing one
#'   combination of these columns' values is handed to \code{fit_fun} as a
#'   single call
#' @param fit_fun A fitting function taking a single group's data frame as
#'   its first argument, returning either a one-row data frame (e.g.
#'   \code{fitRvT_Ea}, \code{fit_temperature_response}) or a named numeric
#'   vector (e.g. \code{fitRvT}, \code{fit_arrhenius}, \code{fitAvT})
#' @param ... Additional arguments passed on to \code{fit_fun}
#'
#' @return A data frame with one row per successfully-fitted group: the
#'   \code{group_cols} values, plus \code{fit_fun}'s result columns. A group
#'   \code{fit_fun} threw an error on (e.g. too few curves for a mixed
#'   model's random effect) is dropped rather than fatally stopping the rest
#'   of the batch, and reported by name with \code{message()} -- the same
#'   \code{tryCatch()} pattern used by this package's own examples of
#'   fitting every group in a dataset, just generalized over which column(s)
#'   define a group
#' @export
fit_Tresponse_groups = function(
  data,
  group_cols,
  fit_fun,
  ...
)
{
  grouped = data |> dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)))
  groups = dplyr::group_keys(grouped)
  data_list = dplyr::group_split(grouped)

  fits_list = lapply(seq_along(data_list), function(i) {
    group_label = paste(paste(group_cols, unlist(groups[i, ])), collapse = ", ")

    tryCatch(
      fit_fun(as.data.frame(data_list[[i]]), ...),
      error = function(e) {
        message("fit_fun() failed for ", group_label, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  ok = !sapply(fits_list, is.null)

  if (!any(ok)) {
    warning("fit_Tresponse_groups(): every group failed to fit")
    return(groups[ok, , drop = FALSE])
  }

  fits_df = do.call(rbind, lapply(fits_list[ok], function(x) {
    if (is.data.frame(x)) x else as.data.frame(t(x))
  }))

  cbind(groups[ok, , drop = FALSE], fits_df, row.names = NULL)
}

# Sensitivity analysis ---------------------------------------------------------

#' Fit a regression to estimate a parameter's climate sensitivity
#'
#' Fits either \code{dpar/dMAT} (a temperature-response parameter's
#' sensitivity to mean annual temperature, when \code{x == "MAT"}) or
#' \code{dToptA/dpar} (assimilation's Topt sensitivity to a perturbed
#' parameter, otherwise), weighted by the parameter's standard error where
#' available.
#'
#' @param x Predictor: \code{"MAT"}, or the name of a perturbed parameter
#'   (matched against \code{df$par})
#' @param y Response column name in \code{df}
#' @param df A data frame of temperature-response fits (when \code{x ==
#'   "MAT"}) or, for a consuming project's own sensitivity-analysis results
#'   (otherwise), a data frame with a \code{par} column matching \code{x}
#'
#' @return A named character/numeric vector: \code{x}, \code{y}, fitted
#'   \code{intercept}, and fitted \code{slope}
#' @export
fit_regression = function(
  x = "MAT",
  y = c("Vcmax25", "EaV", "delsV", "ToptV",
        "Jmax25", "EaJ", "delsJ", "ToptJ"),
  df
)
{
  if (x == "MAT") {
    model = lm(as.formula(paste(y, "~", x)), data = df, weights = 1 / get(paste0(y, ".se"))^2)
  } else {
    df = df |> filter(par == x)
    model = lm(as.formula(paste(y, "~ value")), data = df)
  }
  model.fits = coef(model)
  out = c(x, y, model.fits)
  names(out) = c("x", "y", "intercept", "slope")
  out[3:4] = as.numeric(out[3:4])

  return(out)
}

#' Larger root of a quadratic equation
#'
#' From \pkg{plantecophys} (unpublished function).
#'
#' @param A,B,C Quadratic coefficients (\code{A*x^2 + B*x + C = 0})
#'
#' @return The larger root, or 0 with a warning if the roots are imaginary
#' @export
QUADP <- function(A, B, C)
{
  if ((B^2 - 4*A*C) < 0) {
    warning("IMAGINARY ROOTS IN QUADRATIC")
    return(0)
  }

  if (identical(A, 0)) {
    if (identical(B, 0)) {
      return(0)
    } else {
      return(-C/B)
    }
  } else {
    return((- B + sqrt(B^2 - 4*A*C)) / (2*A))
  }
}

