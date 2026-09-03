# Builds `aci_example`, a small worked-example A-Ci dataset bundled with the
# package (see ?aci_example and the README), from Kumarathunge et al.'s
# published ACi-TGlob database:
#
#   Kumarathunge, D.P. et al. (2019) A-Ci-TGlob: A global dataset of
#   photosynthetic CO2 response curves of terrestrial plants. figshare.
#   https://doi.org/10.6084/m9.figshare.7283567
#
# The full database (~62,000 rows across many species/studies) is too large
# to bundle; this script downloads it fresh and subsets it down to a single
# species/site (Eucalyptus tereticornis, whole-tree chambers, Richmond NSW —
# a Hawkesbury Institute site) that exercises the package's full A-Ci +
# temperature-response workflow: 3 individuals ("replicate"), each with a
# full 5-temperature A-Ci curve series (~20-42C leaf temperature).
#
# Re-run this script (`source("data-raw/aci_example.R")`) to rebuild
# `data/aci_example.rda` if a newer database version changes the subset.

library(dplyr)

url <- "https://ndownloader.figshare.com/files/13449305"
raw_csv <- tempfile(fileext = ".csv")
download.file(url, raw_csv, mode = "wb")

# The source file is not UTF-8 encoded (accented author names in Reference)
aci_tglob <- read.csv(raw_csv, na.strings = "NA", fileEncoding = "latin1")

eter <- aci_tglob |>
  filter(
    Species == "Eucalyptus tereticornis",
    Dataset == "Eucalyptus tereticornis, AU-NSW",
    Temp_Treatment == "ambient"
  )

# Each individual was measured across 5 consecutive Curve_numbers (one per
# target leaf temperature); find runs of 5 consecutive Curve_numbers with
# enough points per curve and a full temperature span
curve_summary <- eter |>
  group_by(Curve_Id, Curve_number) |>
  summarise(n = n(), Tleaf_mean = mean(Tleaf), .groups = "drop") |>
  arrange(Curve_number)

is_run_start <- function(i) {
  five <- curve_summary[i:(i + 4), ]
  all(diff(five$Curve_number) == 1) &&
    all(five$n >= 8) &&
    (max(five$Tleaf_mean) - min(five$Tleaf_mean)) > 15
}
run_starts <- curve_summary$Curve_number[
  vapply(seq_len(nrow(curve_summary) - 4), is_run_start, logical(1))
]

# Three non-overlapping runs, spread across the dataset, as 3 "individuals"
chosen_starts <- c(
  run_starts[1],
  run_starts[run_starts >= run_starts[1] + 20][1],
  run_starts[run_starts >= run_starts[run_starts >= run_starts[1] + 20][1] + 20][1]
)

runs <- lapply(seq_along(chosen_starts), function(i) {
  eter |>
    filter(Curve_number %in% chosen_starts[i]:(chosen_starts[i] + 4)) |>
    mutate(replicate = i)
})

# Label each curve by individual and its rounded target leaf temperature,
# e.g. "rep1_25C"
aci_example <- bind_rows(runs) |>
  group_by(replicate, Curve_number) |>
  mutate(Tleaf_curve_mean = mean(Tleaf)) |>
  ungroup() |>
  mutate(curve_ID = paste0("rep", replicate, "_", round(Tleaf_curve_mean), "C")) |>
  transmute(
    curve_ID,
    replicate,
    species = Species,
    Ci,
    A = Photo,
    Tleaf,
    Qin = PARi,
    Cond
  ) |>
  arrange(replicate, curve_ID, Ci)

usethis::use_data(aci_example, overwrite = TRUE)
