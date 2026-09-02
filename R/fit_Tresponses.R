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

# Fit Vcmax/Jmax temperature responses ----------------------------------------

#' Fit an Arrhenius/peaked-Arrhenius temperature response
#'
#' Fits a peaked Arrhenius model (with an entropy term) and a standard
#' Arrhenius model to \code{Vcmax} or \code{Jmax} vs. leaf temperature, with
#' \code{replicate} as a random effect, and compares the two via AIC/anova.
#' Replaces the previously separate \code{fitpeaked()} (Vcmax) and
#' \code{fitpeakedJ()} (Jmax) functions.
#'
#' @param dat A data frame with columns \code{TsK} (leaf temperature, K),
#'   \code{replicate}, and \code{response}
#' @param response Which response variable to fit: \code{"Vcmax"} or
#'   \code{"Jmax"}
#' @param return Which model's parameters to return: \code{"Arr"} (standard
#'   Arrhenius) or \code{"Peak"} (peaked Arrhenius, plus its fitted Topt)
#' @param start Starting values for the peaked-Arrhenius fit
#'
#' @return A one-row data frame of fitted parameters and standard errors,
#'   plus \code{R2_Arr}, \code{R2_Peak}, \code{AIC_Arr}, \code{AIC_Peak}, and
#'   the anova p-value comparing the two models
#' @export
fit_temperature_response = function(
  dat,
  response = c("Vcmax", "Jmax"),
  return = c("Arr", "Peak"),
  start = c(k25 = 100, Ea = 60, delS = 0.64)
)
{
  response = match.arg(response)
  return = match.arg(return)
  suffix = substr(response, 1, 1)

  dat = dat |> filter(!is.na(!!sym(response)))

  # Peaked Arrhenius model
  Vc <- nlme(model = arrhenius_formula(response, "delS_peaked"),
             fixed = (k25 + Ea + delS ~ 1),
             random = (k25 ~ 1 | replicate),
             data = dat,
             start = start)

  Vc2 <- summary(Vc)
  res1 <- data.frame(t(Vc2$tTable[1:6]))
  names(res1) <- c(
    paste0(response, "25"), paste0("Ea", suffix), paste0("dels", suffix),
    paste0(response, "25.se"), paste0("Ea", suffix, ".se"), paste0("dels", suffix, ".se")
  )

  # Standard Arrhenius model
  Vr <- nlme(model = arrhenius_formula(response, "standard"),
             fixed = (k25 + Ea ~ 1),
             random = (k25 ~ 1 | replicate),
             data = dat,
             start = c(k25 = 100, Ea = 40))

  Vr2 <- summary(Vr)
  res <- data.frame(t(Vr2$tTable[1:4]))
  names(res) <- c(
    paste0(response, "25"), paste0("Ea", suffix),
    paste0(response, "25.se"), paste0("Ea", suffix, ".se")
  )

  an <- anova(Vr, Vc)
  AIC_Arr <- AIC(Vr)
  AIC_Peak <- AIC(Vc)
  prob <- an[[6]][[2]]

  r1 <- cor(fitted(Vc), dat[[response]])
  R2_Peak <- r1 * r1

  r2 <- cor(fitted(Vr), dat[[response]])
  R2_Arr <- r2 * r2

  # Test for normality of residuals (diagnostic only, not returned)
  shapiro.test(residuals(Vc))

  # Topt of the response, with SE (peaked form parameterized by Topt directly)
  VcTopt <- nlme(model = arrhenius_formula(response, "topt_peaked"),
                 fixed = (Ropt + Ea + Topt ~ 1),
                 random = (Ropt ~ 1 | replicate),
                 data = dat,
                 start = c(Ropt = 150, Ea = 60, Topt = 310))

  VcTopt2 <- summary(VcTopt)
  res2 <- data.frame(cbind(VcTopt2$tTable[1], VcTopt2$tTable[3] - 273.15,
                           VcTopt2$tTable[4], VcTopt2$tTable[6]))
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
#' \code{Jmax} vs. leaf temperature, with \code{replicate} as a random
#' effect — simpler than \code{\link{fit_temperature_response}}\code{(return
#' = "Arr")}, which also fits the peaked model for comparison. Replaces the
#' previously separate \code{fitarh()} (Vcmax) and \code{fitarhJ()} (Jmax)
#' functions.
#'
#' @param dat A data frame with columns \code{TsK} (leaf temperature, K),
#'   \code{replicate}, and \code{response}
#' @param response Which response variable to fit: \code{"Vcmax"} or
#'   \code{"Jmax"}
#'
#' @return A named numeric vector of fitted parameters, standard errors, and
#'   \code{R2_Arr}
#' @export
fit_arrhenius = function(
  dat,
  response = c("Vcmax", "Jmax")
)
{
  response = match.arg(response)
  suffix = substr(response, 1, 1)

  dat = dat |> filter(!is.na(!!sym(response)))

  Vr <- nlme(model = arrhenius_formula(response, "standard"),
             fixed = (k25 + Ea ~ 1),
             random = (k25 ~ 1 | replicate),
             data = dat,
             start = c(k25 = 100, Ea = 40))

  Vr2 <- summary(Vr)
  res <- Vr2$tTable[1:4]
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
#' without a per-replicate random effect.
#'
#' @param df A data frame with columns \code{Tleaf}, \code{replicate}, and
#'   \code{var}
#' @param var Which response variable to fit: \code{"A250"}, \code{"Ac250"},
#'   \code{"Aj250"}, \code{"A275"}, or the fully dynamic \code{"A"}
#' @param incl_random_effect Whether to fit a per-\code{replicate} random
#'   effect on \code{Aopt} (via \code{nlme}) or a fixed-effects-only fit
#'   (via \code{nls})
#'
#' @return A named numeric vector of fitted parameters, standard errors, and
#'   the temperature span (\code{Tspan}) for a target fraction of Aopt (see
#'   \code{\link{solve_Tspan}})
#' @export
fitA250vT = function(
  df,
  var = c("A250", "Ac250", "Aj250", "A275", "A"),
  incl_random_effect = TRUE
)
{
  var = match.arg(var)
  df = df |> filter(!is.na(!!sym(var)))

  formula_v = quadratic_topt_formula(var)
  start_vals = c(Aopt = max(df[[var]]), Topt = 25, b = 0.05)

  if (isFALSE(incl_random_effect)) {
    A_Topt <- nls(formula_v, data = df, start = as.list(start_vals))
  } else {
    A_Topt <- nlme(model = formula_v,
                   fixed = (Aopt + b + Topt ~ 1),
                   random = (Aopt ~ 1 | replicate),
                   data = df,
                   start = start_vals)
  }

  A_Topt2 <- summary(A_Topt)

  if (isFALSE(incl_random_effect)) {
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
#'   "MAT"}) or sensitivity-analysis results (otherwise, see
#'   \code{\link{run_sensitivity_analysis}})
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

#' Sensitivity of assimilation's Topt to one Vcmax/Jmax parameter
#'
#' Perturbs a single Vcmax or Jmax temperature-response parameter to a test
#' value, recomputes net assimilation across the curve's measured
#' temperatures, and refits the quadratic Topt model to estimate how much
#' the perturbation shifts the Topt of Rubisco-limited (\code{Ac}),
#' RuBP-regeneration-limited (\code{Aj}), and net (\code{An}) assimilation.
#'
#' @param Tresp_fits A one-row data frame of fitted temperature-response
#'   parameters for one sample (e.g. \code{Vcmax25}, \code{EaV},
#'   \code{delsV}, \code{Jmax25}, \code{EaJ}, \code{delsJ})
#' @param ACi_params A data frame of fitted A-Ci curve parameters for the
#'   same sample (e.g. from \code{\link{extract_curve_params}})
#' @param par Which parameter to perturb: one of \code{"EaV"},
#'   \code{"delsV"}, \code{"Vcmax25"}, \code{"ToptV"}, \code{"EaJ"},
#'   \code{"delsJ"}, \code{"Jmax25"}, \code{"ToptJ"}
#' @param val Test value for \code{par}
#' @param Ci Ci value to evaluate at (a fixed value, or \code{"variable"} to
#'   use each curve's measured Ci — needed for the fully dynamic Topt of A)
#'
#' @return A one-row data frame: \code{par}, \code{value}, and the resulting
#'   \code{Topt_An}/\code{Topt_Ac}/\code{Topt_Aj} (with standard errors)
#' @export
sensitivity_analysis_ToptA = function(
    Tresp_fits,
    ACi_params,
    par = c("EaV", "delsV", "Vcmax25", "ToptV",
            "EaJ", "delsJ", "Jmax", "ToptJ"),
    val,
    Ci = 250
) {
  if (par %in% c("EaV", "delsV", "Vcmax25", "ToptV")) {
    ACi_params = ACi_params |> filter(!(is.na(Jmax) | is.na(Aj250)))
    # Set parameter value
    if (par == "EaV") {
      Ea = val
      k25 = Tresp_fits$Vcmax25
      delS = Tresp_fits$delsV

    } else if (par == "delsV") {
      delS = val
      Ea = Tresp_fits$EaV
      k25 = Tresp_fits$Vcmax25

    } else if (par == "Vcmax25") {
      k25 = val
      Ea = Tresp_fits$EaV
      delS = Tresp_fits$delsV
    } else if (par == "ToptV") {
      Ropt = Tresp_fits$optV
      Ea = Tresp_fits$EaV
      Topt = val + 273.15
    }

    # Fit Vcmax-T response
    TsK = ACi_params$TsK
    Vcmax =
      if (par != "ToptV") { k25 * exp((Ea * (TsK - 298.15))/(298.15 * 0.008314 * TsK)) *
      (1 + exp((298.15 * delS - 200)/(298.15 * 0.008314)))/
      (1 + exp((TsK * delS - 200)/(TsK * 0.008314)))
      } else {
        Ropt*
          ((200*exp(Ea*(TsK-Topt)/(TsK*0.008314*Topt))/
              (200-(Ea*(1-exp(200*(TsK-Topt)/(TsK*0.008314*Topt)))))))
      }

    # Calculate associated net Ac
    if (Ci == "variable") {
      Ci = ACi_params$Ci
    }
    Ac_net = Vcmax * (Ci - ACi_params$GammaStar) /
      (Ci + ACi_params$Km) - ACi_params$Rd

    # Calculate ToptAc
    df = data.frame(A = Ac_net, Tleaf = ACi_params$Tleaf)
    AcvT_fit = fitA250vT(df, var = "A", incl_random_effect = FALSE)
    ToptAc = AcvT_fit[2]
    ToptAc.se = AcvT_fit[5]

    # Calculate net Aj and ToptAj
    Aj_net = ACi_params$Aj250 - ACi_params$Rd
    df = data.frame(A = Aj_net, Tleaf = ACi_params$Tleaf)
    AjvT_fit = fitA250vT(df, var = "A", incl_random_effect = FALSE)
    ToptAj = AjvT_fit[2]
    ToptAj.se = AjvT_fit[5]

  } else if (par %in% c("EaJ", "delsJ", "Jmax25", "ToptJ")) {

    ACi_params = ACi_params |> filter(!(is.na(Vcmax) | is.na(Ac250)))

    # Set parameter value
    if (par == "EaJ") {
      Ea = val
      k25 = Tresp_fits$Jmax25
      delS = Tresp_fits$delsJ

    } else if (par == "delsJ") {
      delS = val
      Ea = Tresp_fits$EaJ
      k25 = Tresp_fits$Jmax25

    } else if (par == "Jmax25") {
      k25 = val
      Ea = Tresp_fits$EaJ
      delS = Tresp_fits$delsJ
    } else if (par == "ToptJ") {
      Ropt = Tresp_fits$optJ
      Ea = Tresp_fits$EaJ
      Topt = val + 273.15
    }

    # Fit Jmax-T response
    TsK = ACi_params$TsK
    Jmax =
      if (par != "ToptJ") {
        k25 * exp((Ea * (TsK - 298.15))/(298.15 * 0.008314 * TsK)) *
          (1 + exp((298.15 * delS - 200)/(298.15 * 0.008314)))/
          (1 + exp((TsK * delS - 200)/(TsK * 0.008314)))
      } else {
        Ropt*
          ((200*exp(Ea*(TsK-Topt)/(TsK*0.008314*Topt))/
              (200-(Ea*(1-exp(200*(TsK-Topt)/(TsK*0.008314*Topt)))))))
      }

    # Calculate associated net Aj
    if (Ci == "variable") {
      Ci = ACi_params$Ci
    }
    Aj_net = Jmax * (Ci - ACi_params$GammaStar) /
      (Ci + 2. * ACi_params$GammaStar) - ACi_params$Rd

    # Calculate ToptAj
    df = data.frame(A = Aj_net, Tleaf = ACi_params$Tleaf)
    AjvT_fit = fitA250vT(df, var = "A", incl_random_effect = FALSE)
    ToptAj = AjvT_fit[2]
    ToptAj.se = AjvT_fit[5]

    # Calculate net Ac and ToptAc
    Ac_net = ACi_params$Ac250 - ACi_params$Rd
    df = data.frame(A = Ac_net, Tleaf = ACi_params$Tleaf)
    AcvT_fit = fitA250vT(df, var = "A", incl_random_effect = FALSE)
    ToptAc = AcvT_fit[2]
    ToptAc.se = AcvT_fit[5]
  }

  Anet <- -mapply(QUADP, A = 1 - 1e-04, B = Ac_net + Aj_net, C = Ac_net *
                  Aj_net)

  # Fit A-T response
  df = data.frame(A = Anet, Tleaf = ACi_params$Tleaf)
  AvT_fit = fitA250vT(df, var = "A", incl_random_effect = FALSE)
  ToptA = AvT_fit[2]
  ToptA.se = AvT_fit[5]

  Topts = data.frame("par" = par, "value" = val,
                     "Topt_An" = ToptA, "Topt_Ac" = ToptAc, "Topt_Aj" = ToptAj,
                     "Topt_An.se" = ToptA.se,"Topt_Ac.se" = ToptAc.se, "Topt_Aj.se" = ToptAj.se
                     )

  return(Topts)
}

#' Run a Vcmax/Jmax parameter sensitivity analysis for one sample
#'
#' For one location/species combination, perturbs each Vcmax/Jmax
#' temperature-response parameter across a small range around its fitted
#' value (see \code{\link{safe_sensitivity_run}}) and collects the resulting
#' shifts in Topt of assimilation.
#'
#' @param all_params_by_sample A named list of fitted A-Ci curve parameter
#'   data frames, one per \code{"<location> - <species>"} combination (e.g.
#'   from splitting \code{\link{extract_curve_params}}'s output)
#' @param Tresp_predictions A data frame of fitted temperature-response
#'   parameters with \code{species} and \code{location} columns
#' @param loc Location name
#' @param sp Species name
#' @param incl_Topts Whether to also perturb \code{ToptV}/\code{ToptJ}
#'   (requires \code{Tresp_predictions} to have those columns)
#'
#' @return A data frame with one row per (parameter, test value) combination
#'   and the resulting \code{Topt_An}/\code{Topt_Ac}/\code{Topt_Aj}
#' @export
run_sensitivity_analysis = function(
    all_params_by_sample,
    Tresp_predictions,
    loc = "Nightcap",
    sp = "Quintinia sieberi",
    incl_Topts = FALSE) {

  # Filter data for location and species of interest
  loc_sp = paste0(loc, " - ", sp)
  ACi_params = all_params_by_sample[[loc_sp]]
  Tresp_fits = Tresp_predictions |>
    filter(species == sp, location == loc)

  # Create data frame with test values of Vcmax/Jmax parameters spanning a
  # small range around the sample's own fitted value
  test_vals = data.frame(EaV = seq(Tresp_fits$EaV - 0.01*Tresp_fits$EaV, Tresp_fits$EaV + 0.01*Tresp_fits$EaV, length.out = 10),
                         EaJ = seq(Tresp_fits$EaJ - 0.01*Tresp_fits$EaJ, Tresp_fits$EaJ + 0.01*Tresp_fits$EaJ, length.out = 10),
                         delsV = seq(Tresp_fits$delsV - 0.01*Tresp_fits$delsV, Tresp_fits$delsV + 0.01*Tresp_fits$delsV, length.out = 10),
                         delsJ = seq(Tresp_fits$delsJ - 0.01*Tresp_fits$delsJ, Tresp_fits$delsJ + 0.01*Tresp_fits$delsJ, length.out = 10),
                         Vcmax25 = seq(Tresp_fits$Vcmax25 - 0.01*Tresp_fits$Vcmax25, Tresp_fits$Vcmax25 + 0.01*Tresp_fits$Vcmax25, length.out = 10),
                         Jmax25 = seq(Tresp_fits$Jmax25 - 0.01*Tresp_fits$Jmax25, Tresp_fits$Jmax25 + 0.01*Tresp_fits$Jmax25, length.out = 10)
  )

  if (isTRUE(incl_Topts)) {
    test_vals$ToptV = c(26.6, Tresp_fits$ToptV, 50.9)
    test_vals$ToptJ = c(19.2, Tresp_fits$ToptJ, 37.8)
  }
  # Loop through each parameter column
  out = lapply(seq_along(test_vals), function(i) {
    par_name = names(test_vals)[i]
    param_vals = test_vals[[i]]

    preds = data.frame(t(sapply(param_vals, function(val) {
      safe_sensitivity_run(par_name, val, Tresp_fits, ACi_params)
    })))

    preds$param = par_name
    preds$test_val = param_vals
    return(preds)
  })

  # Output dataframe reformatting
  out = bind_rows(out)
  out[2:5] = as.numeric(unlist(out[2:5]))
  out[1] = unlist(out[1])

  return(out)
}

#' Solve partial derivatives dA/dMAT from Vcmax/Jmax parameter sensitivities
#'
#' Combines a sample's Vcmax/Jmax parameter sensitivity analysis (see
#' \code{\link{run_sensitivity_analysis}}) with each parameter's fitted
#' climate sensitivity (\code{dpar/dMAT}) to estimate that parameter's
#' contribution to the shift in Topt of assimilation across the climate
#' gradient.
#'
#' @param all_params_by_sample A named list of fitted A-Ci curve parameter
#'   data frames, passed through to \code{\link{run_sensitivity_analysis}}
#' @param Tresp_predictions A data frame of fitted temperature-response
#'   parameters, passed through to \code{\link{run_sensitivity_analysis}}
#' @param par_vs_MAT_fits A data frame of \code{dpar/dMAT} regression fits
#'   (one row per parameter in \code{par_names}, see
#'   \code{\link{fit_regression}}) with a \code{slope} column
#' @param par_names Character vector of Vcmax/Jmax parameter names to
#'   include (e.g. \code{c("EaV", "delsV", "Vcmax25", "EaJ", "delsJ",
#'   "Jmax25")})
#' @param loc Location name
#' @param sp Species name
#' @param incl_Topts Whether to also include \code{ToptV}/\code{ToptJ}
#'
#' @return A data frame with one row per parameter: \code{dpar_dMAT},
#'   \code{dA250_dpar}/\code{dAc_dpar}/\code{dAj_dpar}, and each of their
#'   products with \code{dpar_dMAT}
#' @export
solve_partials = function(
    all_params_by_sample,
    Tresp_predictions,
    par_vs_MAT_fits,
    par_names,
    loc = "Nightcap",
    sp = "Quintinia sieberi",
    incl_Topts = FALSE
    ) {

  # Run sensitivity analysis of Vcmax/Jmax parameters
  out = run_sensitivity_analysis(all_params_by_sample, Tresp_predictions, loc, sp, incl_Topts)

  # Estimate dA/dpar from regression
  ToptA250_vs_par.fits = data.frame(t(sapply(par_names, fit_regression, y = "Topt_An", df = out)))
  ToptA250_vs_par.fits[3:4] = as.numeric(unlist(ToptA250_vs_par.fits[3:4]))

  ToptAc_vs_par.fits = data.frame(t(sapply(par_names, fit_regression, y = "Topt_Ac", df = out)))
  ToptAc_vs_par.fits[3:4] = as.numeric(unlist(ToptAc_vs_par.fits[3:4]))

  ToptAj_vs_par.fits = data.frame(t(sapply(par_names, fit_regression, y = "Topt_Aj", df = out)))
  ToptAj_vs_par.fits[3:4] = as.numeric(unlist(ToptAj_vs_par.fits[3:4]))

  # Create data frame summarizing dA/dpar, dpar/dMAT, and their product
  ToptA250_vs_MAT_props = data.frame(
    location = loc,
    species = sp,
    par = par_names,
    dpar_dMAT = par_vs_MAT_fits$slope,
    dA250_dpar = ToptA250_vs_par.fits$slope,
    A250_partial = ToptA250_vs_par.fits$slope * par_vs_MAT_fits$slope,
    dAc_dpar = ToptAc_vs_par.fits$slope,
    Ac_partial = ToptAc_vs_par.fits$slope * par_vs_MAT_fits$slope,
    dAj_dpar = ToptAj_vs_par.fits$slope,
    Aj_partial = ToptAj_vs_par.fits$slope * par_vs_MAT_fits$slope
  )

  return(ToptA250_vs_MAT_props)
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

#' Run a sensitivity-analysis fit, retrying toward the fitted value on failure
#'
#' Wraps \code{\link{sensitivity_analysis_ToptA}}: if fitting fails to
#' converge at \code{val}, shrinks the gap between \code{val} and the
#' sample's own fitted value by \code{step_fraction} and retries, up to
#' \code{max_attempts} times. This is a deliberate robustness mechanism (test
#' values far from the fitted value more often fail to converge), not error
#' handling to silently bypass failures — if every attempt fails, the result
#' is \code{NA} with a warning.
#'
#' @inheritParams sensitivity_analysis_ToptA
#' @param par_name Which parameter is being perturbed (see \code{par} in
#'   \code{\link{sensitivity_analysis_ToptA}})
#' @param max_attempts Maximum number of retries
#' @param step_fraction Fraction of the remaining gap to \code{Tresp_fits}'s
#'   fitted value to shrink toward on each retry
#'
#' @return The result of \code{\link{sensitivity_analysis_ToptA}}, or a
#'   vector of \code{NA} if every attempt failed
#' @export
safe_sensitivity_run <- function(par_name, val, Tresp_fits, ACi_params, max_attempts = 10, step_fraction = 0.95) {
  # Get the fitted (central) value
  fitted_val <- Tresp_fits[[par_name]]

  attempt <- 1
  success <- FALSE
  result <- NULL

  while (attempt <= max_attempts && !success) {
    try_result <- try({
      out <- sensitivity_analysis_ToptA(val, Tresp_fits = Tresp_fits, ACi_params = ACi_params, par = par_name)
      success <- TRUE
      result <- out
    }, silent = TRUE)

    if (!success) {
      # Shrink the difference between the test value and fitted value
      val <- fitted_val + step_fraction * (val - fitted_val)
      attempt <- attempt + 1
    }
  }

  if (!success) {
    warning(paste("Failed to converge for", par_name, "after", max_attempts, "attempts. Returning NA."))
    return(rep(NA, length(result)))
  }

  return(result)
}

#' @rdname safe_sensitivity_run
#' @export
safe_sensitivity_run_v2 <- function(par_name, val, Tresp_fits, ACi_params, max_attempts = 10, step_fraction = 0.8) {
  fitted_val <- Tresp_fits[[par_name]]

  attempt <- 1
  success <- FALSE
  result <- NULL

  while (attempt <= max_attempts && !success) {
    try_result <- try(
      sensitivity_analysis_ToptA(val, Tresp_fits = Tresp_fits, ACi_params = ACi_params, par = par_name),
      silent = TRUE
    )

    if (!inherits(try_result, "try-error")) {
      success <- TRUE
      result <- try_result
    } else {
      val <- fitted_val + step_fraction * (val - fitted_val)
      attempt <- attempt + 1
    }
  }

  if (!success) {
    warning(paste("Failed to converge for", par_name, "after", max_attempts, "attempts. Returning NA."))
    if (!is.null(result)) {
      return(rep(NA, length(result)))
    } else {
      return(NA)
    }
  }

  return(result)
}
