test_that("add_sp_code derives a 4-letter code from genus + epithet", {
  df = data.frame(species = c("Argyrodendron actinophyllum", "Orites excelsa"))
  out = add_sp_code(df)

  expect_equal(out$sp_code, c("ARAC", "OREX"))
  expect_equal(names(out)[which(names(out) == "species") + 1], "sp_code")
})

test_that("add_species maps codes via a lookup vector, keeping placeholders and NA for unmapped codes", {
  df = data.frame(sp_code = c("ARAC", "BAIN", "ZZZZ"))

  # BAIN is a deliberate placeholder (not yet matched to a taxonomic name),
  # ZZZZ is absent from the lookup entirely
  lookup = c(ARAC = "Argyrodendron actinophyllum", BAIN = "")
  out = add_species(df, lookup)

  expect_equal(out$species, c("Argyrodendron actinophyllum", "", NA_character_))
})

test_that("add_species also accepts a 2-column data frame lookup", {
  df = data.frame(sp_code = c("ARAC", "OREX"))
  lookup = data.frame(
    sp_code = c("ARAC", "OREX"),
    species = c("Argyrodendron actinophyllum", "Orites excelsa")
  )

  out = add_species(df, lookup)
  expect_equal(out$species, c("Argyrodendron actinophyllum", "Orites excelsa"))
})

# The fixture file lacks embedded factory calibration data (a characteristic
# of this particular raw export, not of tidy_gasxchg()), which gasanalyzer
# warns about -- expected and harmless for what this test checks (filename
# parsing, column presence), so it's asserted explicitly rather than left as
# a stray warning.
test_that("tidy_gasxchg parses A-Ci filename metadata correctly", {
  withr::local_dir(test_path("fixtures"))

  expect_warning(
    df <- tidy_gasxchg(
      "data/in/raw/20260504/dipe_700m_rep1_aci_25C.xlsx",
      type = "ACi",
      area_correction = NA
    ),
    "factory calibration"
  )

  expect_true(nrow(df) > 0)
  expect_equal(unique(df$curve_ID), "dipe_700m_rep1_aci_25C")
  expect_equal(unique(df$sp_code), "DIPE")
  expect_equal(unique(df$elevation), "700m")
  expect_equal(unique(df$replicate), "rep1")
  expect_equal(unique(df$meas_temp), "25C")
  expect_equal(unique(df$type), "ACi-T")
  expect_true(all(c("A", "Ci", "Tleaf", "Patm", "RHs", "area", "flow") %in% names(df)))
})

test_that("tidy_gasxchg negates A into Rdark for R-T curves", {
  withr::local_dir(test_path("fixtures"))

  df = tidy_gasxchg(
    "data/in/raw/20260504/bain_700m_rep2_rt.xlsx",
    type = "RT",
    area_correction = NA
  )

  expect_true(nrow(df) > 0)
  expect_equal(unique(df$curve_ID), "bain_700m_rep2_rt")
  expect_equal(unique(df$sp_code), "BAIN")
  expect_equal(unique(df$replicate), "rep2")
  expect_equal(unique(df$type), "RT")
  expect_false("A" %in% names(df))
  expect_true("Rdark" %in% names(df))
})
