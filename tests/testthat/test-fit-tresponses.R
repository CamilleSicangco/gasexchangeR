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

# Synthetic temperature-response data with a per-replicate random effect on
# k25, for fit_arrhenius()/fit_temperature_response().
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
  replicate = rep(paste0("R", seq_len(n_rep)), times = length(Tleaf_C))
  rep_effect = stats::setNames(rnorm(n_rep, 0, rep_sd), paste0("R", seq_len(n_rep)))
  k25_i = k25 + rep_effect[replicate]

  y = if (is.null(delS)) {
    k25_i * exp((Ea * (TsK - 298.15)) / (298.15 * 0.008314 * TsK))
  } else {
    k25_i * exp((Ea * (TsK - 298.15)) / (298.15 * 0.008314 * TsK)) *
      (1 + exp((298.15 * delS - 200) / (298.15 * 0.008314))) /
      (1 + exp((TsK * delS - 200) / (TsK * 0.008314)))
  }
  y = y + rnorm(length(y), 0, noise_sd)

  df = data.frame(TsK = TsK, replicate = replicate, y = y)
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

test_that("fitA250vT recovers a known quadratic Topt response", {
  set.seed(24)
  n_rep = 4
  Tleaf = rep(seq(15, 40, by = 5), each = n_rep)
  replicate = rep(paste0("R", seq_len(n_rep)), times = length(unique(Tleaf)))
  Aopt_true = 20
  Topt_true = 28
  b_true = 0.03
  A250 = Aopt_true - b_true * (Tleaf - Topt_true)^2 + rnorm(length(Tleaf), 0, 0.3)
  df = data.frame(Tleaf = Tleaf, replicate = replicate, A250 = A250)

  fit = fitA250vT(df, var = "A250", incl_random_effect = TRUE)

  expect_equal(unname(fit["A250opt"]), Aopt_true, tolerance = 0.1)
  expect_equal(unname(fit["ToptA250"]), Topt_true, tolerance = 0.15)

  fit_fixed = fitA250vT(df, var = "A250", incl_random_effect = FALSE)
  expect_equal(unname(fit_fixed["A250opt"]), Aopt_true, tolerance = 0.1)
})
