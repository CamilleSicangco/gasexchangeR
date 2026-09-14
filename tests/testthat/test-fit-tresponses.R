test_that("arrhenius_formula builds the correct formula for each response/form", {
  # Compared as text: the formulas' *content* is what matters here, not
  # which environment each happens to capture (arrhenius_formula()'s own
  # local variables vs. this test's), which expect_equal() on formula
  # objects would otherwise flag as a difference.
  expect_equal(
    deparse(arrhenius_formula("Vcmax", "standard")),
    deparse(as.formula("Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK))"))
  )
  expect_equal(
    deparse(arrhenius_formula("Jmax", "delS_peaked")),
    deparse(as.formula(paste0(
      "Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * ",
      "(1+exp((298.15*delS - 200)/(298.15*0.008314))) / ",
      "(1+exp((TsK*delS-200)/(TsK*0.008314)))"
    )))
  )
  expect_equal(
    deparse(arrhenius_formula("Vcmax", "topt_peaked")),
    deparse(as.formula(paste0(
      "Vcmax ~ Ropt*((200*exp(Ea*(TsK-Topt)/(TsK*0.008314*Topt))/",
      "(200-(Ea*(1-exp(200*(TsK-Topt)/(TsK*0.008314*Topt)))))))"
    )))
  )
})

test_that("quadratic_topt_formula builds the correct formula for each response", {
  expect_equal(
    deparse(quadratic_topt_formula("A250")),
    deparse(as.formula("A250 ~ Aopt-(b*(Tleaf-Topt)^2)"))
  )
  expect_equal(
    deparse(quadratic_topt_formula("A")),
    deparse(as.formula("A ~ Aopt-(b*(Tleaf-Topt)^2)"))
  )
})

test_that("solve_Tspan matches a hand-computed value", {
  # rt = sqrt((1-0.8)*20/0.05) = sqrt(80) = 8.944272; Tspan = 2*rt
  expect_equal(solve_Tspan(Aopt = 20, Topt = 30, b = 0.05), 2 * sqrt(80))
})

test_that("QUADP returns the larger root of a quadratic, and handles degenerate cases", {
  expect_equal(QUADP(1, -3, 2), 2)   # roots 1, 2
  expect_equal(QUADP(1, -5, 6), 3)   # roots 2, 3
  expect_equal(QUADP(0, 5, 10), -2)  # linear: 5x + 10 = 0
  expect_equal(QUADP(0, 0, 5), 0)
  expect_warning(expect_equal(QUADP(1, 0, 1), 0), "IMAGINARY ROOTS")
})

# Synthetic temperature-response data with a per-curve_ID random effect on
# k25, for fit_arrhenius()/fit_temperature_response() (curve_ID being their
# default random_effect_col).
simulate_arrhenius_data = function(
  k25 = 60,
  Ea = 50,
  delS = NULL,
  n_rep = 4,
  Tleaf_C = c(15, 20, 25, 30, 35, 40, 45),
  noise_sd = 1,
  rep_sd = 2,
  seed = 1,
  response = "Vcmax"
) {
  set.seed(seed)
  TsK = rep(Tleaf_C, each = n_rep) + 273.15
  curve_ID = rep(paste0("R", seq_len(n_rep)), times = length(Tleaf_C))
  rep_effect = stats::setNames(rnorm(n_rep, 0, rep_sd), paste0("R", seq_len(n_rep)))
  k25_i = k25 + rep_effect[curve_ID]

  y = if (is.null(delS)) {
    k25_i * exp((Ea * (TsK - 298.15)) / (298.15 * 0.008314 * TsK))
  } else {
    k25_i * exp((Ea * (TsK - 298.15)) / (298.15 * 0.008314 * TsK)) *
      (1 + exp((298.15 * delS - 200) / (298.15 * 0.008314))) /
      (1 + exp((TsK * delS - 200) / (TsK * 0.008314)))
  }
  y = y + rnorm(length(y), 0, noise_sd)

  df = data.frame(TsK = TsK, curve_ID = curve_ID, y = y)
  names(df)[3] = response
  df
}

test_that("fit_arrhenius recovers known Vcmax25/Ea from synthetic data", {
  dat = simulate_arrhenius_data(k25 = 60, Ea = 50, response = "Vcmax", seed = 21)

  fit = fit_arrhenius(dat, response = "Vcmax")

  expect_equal(unname(fit["Vcmax25"]), 60, tolerance = 0.1)
  expect_equal(unname(fit["EaV"]), 50, tolerance = 0.15)
  expect_true(unname(fit["R2_Arr"]) > 0.9)
})

test_that("fit_arrhenius recovers known Jmax25/Ea from synthetic data", {
  dat = simulate_arrhenius_data(k25 = 120, Ea = 40, response = "Jmax", seed = 22)

  fit = fit_arrhenius(dat, response = "Jmax")

  expect_equal(unname(fit["Jmax25"]), 120, tolerance = 0.1)
  expect_equal(unname(fit["EaJ"]), 40, tolerance = 0.15)
})

test_that("fit_arrhenius fits without error using a custom random_effect_col", {
  dat = simulate_arrhenius_data(k25 = 60, Ea = 50, response = "Vcmax", seed = 21)
  names(dat)[names(dat) == "curve_ID"] = "individual"

  fit = fit_arrhenius(dat, response = "Vcmax", random_effect_col = "individual")

  expect_equal(unname(fit["Vcmax25"]), 60, tolerance = 0.1)
})

test_that("fit_arrhenius falls back to a fixed-effects fit when random_effect_col has one level", {
  dat = simulate_arrhenius_data(k25 = 60, Ea = 50, response = "Vcmax", seed = 21)
  dat$curve_ID = "only_one_curve"

  expect_message(
    fit <- fit_arrhenius(dat, response = "Vcmax"),
    "Only one level of 'curve_ID'"
  )

  expect_false(anyNA(fit))
  expect_true(all(c("Vcmax25", "EaV", "R2_Arr") %in% names(fit)))
})

test_that("fit_temperature_response recovers known peaked-Arrhenius parameters", {
  dat = simulate_arrhenius_data(k25 = 60, Ea = 50, delS = 0.64, response = "Vcmax",
                                 noise_sd = 0.5, seed = 23)

  fit = fit_temperature_response(dat, response = "Vcmax", return = "Peak")

  expect_equal(unname(fit[["Vcmax25"]]), 60, tolerance = 0.15)
  expect_equal(unname(fit[["EaV"]]), 50, tolerance = 0.2)
  expect_true(unname(fit[["R2_Peak"]]) > 0.9)

  fit_arr = fit_temperature_response(dat, response = "Vcmax", return = "Arr")
  expect_true(all(c("Vcmax25", "EaV", "R2_Arr", "R2_Peak", "AIC_Arr", "AIC_Peak", "prob_V") %in% names(fit_arr)))
})

test_that("fit_temperature_response falls back to fixed-effects fits when random_effect_col has one level", {
  dat = simulate_arrhenius_data(k25 = 60, Ea = 50, delS = 0.64, response = "Vcmax",
                                 noise_sd = 0.5, seed = 23)
  dat$curve_ID = "only_one_curve"

  expect_message(
    fit <- fit_temperature_response(dat, response = "Vcmax", return = "Peak"),
    "Only one level of 'curve_ID'"
  )

  expect_false(anyNA(fit))
  expect_true(all(c("Vcmax25", "EaV", "R2_Arr", "R2_Peak", "AIC_Arr", "AIC_Peak", "prob_V") %in% names(fit)))
})

test_that("fitRvT recovers a known Q10-type respiration-temperature response", {
  set.seed(25)
  n_rep = 4
  Tleaf = rep(seq(10, 35, by = 5), each = n_rep)
  sample = rep(paste0("S", seq_len(n_rep)), times = length(unique(Tleaf)))
  a_true = log(2); b_true = 0.05; c_true = -0.0005
  rep_effect = stats::setNames(rnorm(n_rep, 0, 0.05), paste0("S", seq_len(n_rep)))
  ln_R = a_true + rep_effect[sample] + b_true * Tleaf + c_true * Tleaf^2 +
    rnorm(length(Tleaf), 0, 0.03)
  df = data.frame(Tleaf = Tleaf, sample = sample, R = exp(ln_R))

  fit = fitRvT(df, random_effect_col = "sample")

  Rdark25_true = exp(a_true + b_true * 25 + c_true * 25^2)
  expect_equal(unname(fit["Rdark25"]), Rdark25_true, tolerance = 0.15)
  expect_true(unname(fit["r2"]) > 0.9)

  # Same data, differently-named respiration column, selected via var=
  names(df)[names(df) == "R"] = "Rday"
  fit_rday = fitRvT(df, var = "Rday", random_effect_col = "sample")
  expect_equal(unname(fit_rday["Rdark25"]), Rdark25_true, tolerance = 0.15)

  # Same data again, differently-named grouping column, selected via
  # random_effect_col=
  names(df)[names(df) == "sample"] = "individual"
  fit_custom_col = fitRvT(df, var = "Rday", random_effect_col = "individual")
  expect_equal(unname(fit_custom_col["Rdark25"]), Rdark25_true, tolerance = 0.15)
})

test_that("fitRvT falls back to a fixed-effects fit when random_effect_col has one level", {
  set.seed(26)
  Tleaf = seq(10, 35, by = 5)
  a_true = log(2); b_true = 0.05; c_true = -0.0005
  ln_R = a_true + b_true * Tleaf + c_true * Tleaf^2 + rnorm(length(Tleaf), 0, 0.03)
  df = data.frame(Tleaf = Tleaf, curve_ID = "only_one_curve", R = exp(ln_R))

  expect_message(
    fit <- fitRvT(df),
    "Only one level of 'curve_ID'"
  )

  Rdark25_true = exp(a_true + b_true * 25 + c_true * 25^2)
  expect_equal(unname(fit["Rdark25"]), Rdark25_true, tolerance = 0.15)
})

test_that("fit_Tresponse_groups fits every group, attaches group columns, and drops failures", {
  # Mirrors the fitRvT test's own data-generating pattern above, just
  # parameterized by n_curves (its n_rep) so a single-curve group can be
  # simulated too.
  simulate_one_group = function(a_true, b_true, c_true, n_curves, seed) {
    set.seed(seed)
    Tleaf_seq = seq(10, 35, by = 5)
    Tleaf_g = rep(Tleaf_seq, each = n_curves)
    sample = rep(paste0("S", seq_len(n_curves)), times = length(Tleaf_seq))
    rep_effect = stats::setNames(rnorm(n_curves, 0, 0.05), paste0("S", seq_len(n_curves)))
    ln_R = a_true + rep_effect[sample] + b_true * Tleaf_g + c_true * Tleaf_g^2 +
      rnorm(length(Tleaf_g), 0, 0.03)
    data.frame(Tleaf = Tleaf_g, sample = sample, R = exp(ln_R))
  }

  # Two groups with enough curves for fitRvT's per-sample random effect, one
  # with only a single curve (no replication) -- exactly like a real
  # single-curve species/elevation group. All three should still fit, the
  # last one falling back to a fixed-effects-only fit.
  group_a = cbind(sp_code = "AAAA", elevation = "500m",
                  simulate_one_group(log(2), 0.05, -0.0005, n_curves = 4, seed = 261))
  group_b = cbind(sp_code = "BBBB", elevation = "500m",
                  simulate_one_group(log(3), 0.04, -0.0003, n_curves = 4, seed = 262))
  group_c = cbind(sp_code = "AAAA", elevation = "900m",
                  simulate_one_group(log(1.5), 0.06, -0.0004, n_curves = 1, seed = 263))

  all_data = rbind(group_a, group_b, group_c)

  expect_message(
    result <- fit_Tresponse_groups(all_data, group_cols = c("sp_code", "elevation"),
                                    fit_fun = fitRvT, random_effect_col = "sample"),
    "Only one level of 'sample'"
  )

  expect_equal(nrow(result), 3)
  expect_true(all(c("sp_code", "elevation", "Rdark25", "r2") %in% names(result)))

  # True Rdark25 = exp(a + b*25 + c*25^2), matching the fitRvT test above
  fit_a = result[result$sp_code == "AAAA" & result$elevation == "500m", ]
  expect_equal(fit_a$Rdark25, exp(log(2) + 0.05 * 25 - 0.0005 * 25^2), tolerance = 0.2)

  fit_b = result[result$sp_code == "BBBB" & result$elevation == "500m", ]
  expect_equal(fit_b$Rdark25, exp(log(3) + 0.04 * 25 - 0.0003 * 25^2), tolerance = 0.2)

  # The single-curve group: fit without a random effect, but still fit
  fit_c = result[result$sp_code == "AAAA" & result$elevation == "900m", ]
  expect_equal(nrow(fit_c), 1)
  expect_equal(fit_c$Rdark25, exp(log(1.5) + 0.06 * 25 - 0.0004 * 25^2), tolerance = 0.2)
})

test_that("fitAvT recovers a known quadratic Topt response", {
  set.seed(24)
  n_rep = 4
  Tleaf = rep(seq(15, 40, by = 5), each = n_rep)
  curve_ID = rep(paste0("R", seq_len(n_rep)), times = length(unique(Tleaf)))
  Aopt_true = 20
  Topt_true = 28
  b_true = 0.03
  A250 = Aopt_true - b_true * (Tleaf - Topt_true)^2 + rnorm(length(Tleaf), 0, 0.3)
  df = data.frame(Tleaf = Tleaf, curve_ID = curve_ID, A250 = A250)

  fit = fitAvT(df, var = "A250", incl_random_effect = TRUE)

  expect_equal(unname(fit["A250opt"]), Aopt_true, tolerance = 0.1)
  expect_equal(unname(fit["ToptA250"]), Topt_true, tolerance = 0.15)

  fit_fixed = fitAvT(df, var = "A250", incl_random_effect = FALSE)
  expect_equal(unname(fit_fixed["A250opt"]), Aopt_true, tolerance = 0.1)

  # Custom random_effect_col name
  names(df)[names(df) == "curve_ID"] = "individual"
  fit_custom_col = fitAvT(df, var = "A250", random_effect_col = "individual")
  expect_equal(unname(fit_custom_col["A250opt"]), Aopt_true, tolerance = 0.1)
})

test_that("fitAvT falls back to a fixed-effects fit when random_effect_col has one level", {
  set.seed(24)
  n_rep = 4
  Tleaf = rep(seq(15, 40, by = 5), each = n_rep)
  Aopt_true = 20
  Topt_true = 28
  b_true = 0.03
  A250 = Aopt_true - b_true * (Tleaf - Topt_true)^2 + rnorm(length(Tleaf), 0, 0.3)
  df = data.frame(Tleaf = Tleaf, curve_ID = "only_one_curve", A250 = A250)

  expect_message(
    fit <- fitAvT(df, var = "A250", incl_random_effect = TRUE),
    "Only one level of 'curve_ID'"
  )

  expect_equal(unname(fit["A250opt"]), Aopt_true, tolerance = 0.1)
})

test_that("fitAvT accepts an arbitrary var/response name outside the known presets", {
  set.seed(27)
  n_rep = 4
  Tleaf = rep(seq(15, 40, by = 5), each = n_rep)
  curve_ID = rep(paste0("R", seq_len(n_rep)), times = length(unique(Tleaf)))
  Aopt_true = 20
  Topt_true = 28
  b_true = 0.03
  A_custom = Aopt_true - b_true * (Tleaf - Topt_true)^2 + rnorm(length(Tleaf), 0, 0.3)
  df = data.frame(Tleaf = Tleaf, curve_ID = curve_ID, A_custom = A_custom)

  fit = fitAvT(df, var = "A_custom", incl_random_effect = TRUE)

  expect_true(all(c("A_customopt", "ToptA_custom", "b_A_custom", "Tspan_A_custom") %in% names(fit)))
  expect_equal(unname(fit["A_customopt"]), Aopt_true, tolerance = 0.1)
  expect_equal(unname(fit["ToptA_custom"]), Topt_true, tolerance = 0.2)
})

test_that("fitAvT's Topt_start defaults to the observed max and accepts a fixed override", {
  set.seed(28)
  n_rep = 4
  Tleaf = rep(seq(15, 40, by = 5), each = n_rep)
  curve_ID = rep(paste0("R", seq_len(n_rep)), times = length(unique(Tleaf)))
  Aopt_true = 20
  Topt_true = 32
  b_true = 0.02
  A250 = Aopt_true - b_true * (Tleaf - Topt_true)^2 + rnorm(length(Tleaf), 0, 0.3)
  df = data.frame(Tleaf = Tleaf, curve_ID = curve_ID, A250 = A250)

  fit_default = fitAvT(df, var = "A250", incl_random_effect = TRUE)
  expect_equal(unname(fit_default["ToptA250"]), Topt_true, tolerance = 0.15)

  fit_fixed_start = fitAvT(df, var = "A250", incl_random_effect = TRUE, Topt_start = 25)
  expect_equal(unname(fit_fixed_start["ToptA250"]), Topt_true, tolerance = 0.15)
})
