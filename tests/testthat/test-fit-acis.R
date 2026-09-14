# Simulate a simple two-limitation (Rubisco- and RuBP-regeneration-limited)
# A-Ci curve from known Vcmax/Jmax/Rd, assuming light-saturated conditions
# (so J ~= Jmax) -- enough to exercise curve fitting without needing a full
# photosynthesis simulator.
simulate_aci_curve = function(
  Ci,
  Vcmax = 80,
  Jmax = 150,
  Rd = 1.2,
  Tleaf = 25,
  GammaStar = 40,
  Km = 700,
  noise_sd = 0.3,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  Ac = Vcmax * (Ci - GammaStar) / (Ci + Km) - Rd
  Aj = Jmax * (Ci - GammaStar) / (4 * (Ci + 2 * GammaStar)) - Rd
  A = pmin(Ac, Aj) + rnorm(length(Ci), 0, noise_sd)
  data.frame(Ci = Ci, A = A, Tleaf = Tleaf, Qin = 1500)
}

test_that("fit_ACis/get_r2/extract_curve_params recover known parameters (plantecophys engine)", {
  Ci_vals = c(50, 100, 150, 200, 300, 400, 600, 800, 1200, 1500)

  curve1 = simulate_aci_curve(Ci_vals, Vcmax = 80, Jmax = 150, Rd = 1.2, seed = 10)
  curve1$curve_ID = "C1"
  curve2 = simulate_aci_curve(Ci_vals, Vcmax = 60, Jmax = 120, Rd = 1.0, seed = 11)
  curve2$curve_ID = "C2"
  dat = rbind(curve1, curve2)

  fits = fit_ACis(dat, group_cols = "curve_ID", engine = "plantecophys")
  expect_length(fits, 2)

  r2 = get_r2(fits, engine = "plantecophys")
  expect_true(all(r2 > 0.9))

  extra = data.frame(
    curve_ID = c("C1", "C2"),
    location = "TestLoc",
    sample = c("S1", "S2"),
    method = "test",
    Tleaf = c(25, 25)
  )
  params = extract_curve_params(fits, extra, byname = "curve_ID", engine = "plantecophys")

  expect_true(all(c("Vcmax", "Jmax", "Rd") %in% names(params)))
  vc = params$Vcmax[match(c("C1", "C2"), params$curve_ID)]
  jm = params$Jmax[match(c("C1", "C2"), params$curve_ID)]
  expect_equal(vc, c(80, 60), tolerance = 0.15)
  expect_equal(jm, c(150, 120), tolerance = 0.15)
})

# The PhotoGEA-engine tests below exercise get_r2()/extract_curve_params()'s
# PhotoGEA branch directly against hand-built fixtures shaped like a real
# PhotoGEA fit result, rather than running fit_ACis(engine = "PhotoGEA") end
# to end. A real PhotoGEA fit requires a full `exdf` (with unit-typed
# columns) and by default uses a slow global optimizer (DEoptim) -- better
# suited to a slower, opt-in integration test than this fast unit-test
# suite. The dispatch/extraction logic these tests cover is what this
# package's engine-unification actually changed; PhotoGEA's own fit_c3_aci()
# correctness is that package's responsibility, not re-tested here.

test_that("get_r2 computes R2 from a PhotoGEA-shaped fit result", {
  fits_main = data.frame(
    curve_ID = c("C1", "C1", "C1", "C2", "C2", "C2"),
    A = c(10, 15, 18, 8, 12, 14),
    A_fit = c(10.2, 14.8, 18.1, 7.9, 12.3, 13.8)
  )
  fit_result = list(fits = list(main_data = fits_main))

  r2 = get_r2(fit_result, engine = "PhotoGEA")

  expect_equal(names(r2), c("C1", "C2"))
  expect_true(all(r2 > 0.95))
})

test_that("extract_curve_params extracts fixed-Ci values and merges data (PhotoGEA engine)", {
  params_main = data.frame(
    curve_ID = c("C1", "C2"),
    TleafCnd_avg = c(25, 26),
    Vcmax_tl_avg = c(80, 60),
    Jmax_tl = c(150, 120),
    RL_tl_avg = c(1.2, 1.0),
    Tp_tl_avg = c(10, 8),
    Gamma_star_tl_avg = c(40, 41)
  )

  fits_main = data.frame(
    curve_ID = c("C1", "C1", "C1", "C2", "C2", "C2"),
    A = c(10, 15, 18, 8, 12, 14),
    A_fit = c(10.2, 14.8, 18.1, 7.9, 12.3, 13.8)
  )

  interp = data.frame(
    curve_ID = rep(c("C1", "C2"), each = 5),
    Ci = rep(c(100, 250, 275, 360, 500), times = 2),
    An = c(10, 18, 19, 20, 15, 8, 14, 15, 16, 12),
    Ac = c(9, 17, 18, 19, 14, 7, 13, 14, 15, 11),
    Aj = c(11, 19, 20, 21, 16, 9, 15, 16, 17, 13),
    Ap = c(12, 20, 21, 22, 17, 10, 16, 17, 18, 14)
  )

  fits_exdf = list(
    parameters = list(main_data = params_main),
    fits = list(main_data = fits_main),
    fits_interpolated = list(main_data = interp)
  )

  data = data.frame(
    curve_ID = c("C1", "C2"),
    species = c("Sp1", "Sp2"),
    E = c(0.002, 0.003),
    Tirga = c(30, 31)
  )

  ret = extract_curve_params(fits_exdf, data, byname = "curve_ID", engine = "PhotoGEA")

  expect_equal(nrow(ret), 2)
  expect_true(all(c("Vcmax", "Jmax", "Rd", "TPU", "GammaStar", "A250", "A275", "A360") %in% names(ret)))

  ret = ret[match(c("C1", "C2"), ret$curve_ID), ]
  expect_equal(ret$Vcmax, c(80, 60))
  expect_equal(ret$A250, c(18, 14))
  expect_equal(ret$species, c("Sp1", "Sp2"))
})

test_that("get_r2 accepts a custom identifier_column name (PhotoGEA engine)", {
  fits_main = data.frame(
    curve_id = c("C1", "C1", "C1", "C2", "C2", "C2"),
    A = c(10, 15, 18, 8, 12, 14),
    A_fit = c(10.2, 14.8, 18.1, 7.9, 12.3, 13.8)
  )
  fit_result = list(fits = list(main_data = fits_main))

  r2 = get_r2(fit_result, engine = "PhotoGEA", identifier_column = "curve_id")

  expect_equal(names(r2), c("C1", "C2"))
  expect_true(all(r2 > 0.95))
})

test_that("extract_curve_params accepts Ci_targets and a differently-named identifier column (PhotoGEA engine)", {
  # Same shape as the "extracts fixed-Ci values and merges data" test above,
  # but with the identifier column named "curve_id" (lowercase, matching a
  # real consumer's own naming) instead of the default "curve_ID", and an
  # extra Ci_targets entry (300) beyond the historical 250/275/360.
  params_main = data.frame(
    curve_id = c("C1", "C2"),
    TleafCnd_avg = c(25, 26),
    Vcmax_tl_avg = c(80, 60),
    Jmax_tl = c(150, 120),
    RL_tl_avg = c(1.2, 1.0),
    Tp_tl_avg = c(10, 8),
    Gamma_star_tl_avg = c(40, 41)
  )

  fits_main = data.frame(
    curve_id = c("C1", "C1", "C1", "C2", "C2", "C2"),
    A = c(10, 15, 18, 8, 12, 14),
    A_fit = c(10.2, 14.8, 18.1, 7.9, 12.3, 13.8)
  )

  interp = data.frame(
    curve_id = rep(c("C1", "C2"), each = 6),
    Ci = rep(c(100, 250, 275, 300, 360, 500), times = 2),
    An = c(10, 18, 19, 19.5, 20, 15, 8, 14, 15, 15.5, 16, 12),
    Ac = c(9, 17, 18, 18.5, 19, 14, 7, 13, 14, 14.5, 15, 11),
    Aj = c(11, 19, 20, 20.5, 21, 16, 9, 15, 16, 16.5, 17, 13),
    Ap = c(12, 20, 21, 21.5, 22, 17, 10, 16, 17, 17.5, 18, 14)
  )

  fits_exdf = list(
    parameters = list(main_data = params_main),
    fits = list(main_data = fits_main),
    fits_interpolated = list(main_data = interp)
  )

  data = data.frame(
    curve_id = c("C1", "C2"),
    species = c("Sp1", "Sp2"),
    E = c(0.002, 0.003),
    Tirga = c(30, 31)
  )

  ret = extract_curve_params(fits_exdf, data, byname = "curve_id", engine = "PhotoGEA",
                             Ci_targets = c(250, 275, 300, 360))

  expect_equal(nrow(ret), 2)
  expect_true(all(c("A250", "A275", "A300", "A360") %in% names(ret)))

  ret = ret[match(c("C1", "C2"), ret$curve_id), ]
  expect_equal(ret$A300, c(19.5, 15.5))
})
