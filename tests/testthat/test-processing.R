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
