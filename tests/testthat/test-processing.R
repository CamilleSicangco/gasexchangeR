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

# The tests below require reading real LI-6800 .xlsx fixtures via
# gasanalyzer::read_6800_xlsx(), which in turn depends on tidyxl. As of
# writing, xlsx reading crashes R entirely in the development environment
# (reproduced with tidyxl::xlsx_cells() and readxl::read_excel() on a
# trivial, freshly-generated .xlsx file, under both R 4.5.2 and R 4.3.1) --
# an apparent system-level issue unrelated to this package. These tests are
# believed correct but could not be executed/verified here; please confirm
# they pass in an environment where xlsx reading works.
test_that("tidy_gasxchg parses A-Ci filename metadata correctly", {
  skip("xlsx reading crashes R in the current dev environment (see comment above) -- remove this skip() once verified to work")
  withr::local_dir(test_path("fixtures"))

  df = tidy_gasxchg(
    "data/in/raw/20260504/dipe_700m_rep1_aci_25C.xlsx",
    type = "ACi",
    area_correction = NA
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
  skip("xlsx reading crashes R in the current dev environment (see comment above) -- remove this skip() once verified to work")
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
