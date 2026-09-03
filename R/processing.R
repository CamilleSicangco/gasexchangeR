#' Fixed column order for tidied LI-6800 output
#'
#' Internal constant listing the columns kept and their order in the output
#' of \code{\link{tidy_gasxchg}}, using ESS-DIVE-style names.
#' @noRd
cols_to_select = c(
  "Date", "HHMMSS", "curve_ID", "sp_code", "elevation", "replicate", "type",
  "meas_temp","Instrument", "UseDynamic", "area", "K", "Oxygen",
  "E", "A", "Rdark", "Ca", "Ci", "gsw", "gbw", "gtc", "gtw", "Rabs",
  "TleafEB", "TleafCnd", "SVPleaf", "RHs", "VPcham", "SVPcham",
  "VPDleaf", "LatHFlux", "SenHFlux", "NetTherm", "EBSum", "Asty",
  "Esty", "Adyn", "Crd", "Csd", "dCsd_dt", "alphaVc", "Edyn", "Hr", "Hs",
  "dHs_dt", "alphaVh", "Qabs", "Qin", "fan", "flow", "Patm", "DeltaPcham",
  "Tair", "Tleaf", "Tleaf2", "FanSpeed", "QambIn", "QambOut", "CO2s",
  "CO2a", "CO2r", "H2Os", "H2Oa", "H2Or"
)

#' Re-apply a manually measured leaf area to a LI-6800 file
#'
#' Reads a raw LI-6800 Excel export, overwrites its embedded leaf-area
#' constant with a manually measured value, and recalculates every derived
#' gas-exchange quantity using that file's own embedded equations.
#'
#' @param file_path Path to a raw LI-6800 \code{.xlsx} export
#' @param corrected_area Manually measured leaf area, cm^2
#'
#' @return An \code{exdf} object (see \code{gasanalyzer::read_6800_xlsx}) with
#'   leaf area and all area-dependent quantities recalculated
#' @export
correct_6800_area = function(
  file_path,
  corrected_area
)
{
  # Read in data
  xlsxdat <- read_6800_xlsx(file_path)

  # A file with no data rows can't be recalculated (transform() errors on a
  # 0-row data frame) - leave it as-is and let the caller (tidy_gasxchg)
  # detect and drop it
  if (nrow(xlsxdat) == 0) {
    return(xlsxdat)
  }

  # Correct the area
  xlsxdat_corrected <- xlsxdat |>
    transform(Const.S = set_units(corrected_area, "cm^2")) |>
    recalculate()   # uses the equations embedded in this specific file

  return(xlsxdat_corrected)
}

#' Tidy a single LI-6800 output file
#'
#' Reads a raw LI-6800 Excel export (optionally re-applying a manually
#' measured leaf area), drops duplicated/raw columns, strips column-group
#' prefixes, parses identifying information out of the filename, renames
#' columns to ESS-DIVE-style names, and selects/orders a fixed set of output
#' columns.
#'
#' The filename (see \code{path}) is expected to follow the convention
#' \verb{<sp_code>_<elevation>_<rep>[_newleaf]_<aci|rt>_<temp>C} for A-Ci
#' curves, or \verb{<sp_code>_<elevation>_<rep>_rt} for R-T curves, with the
#' curve identifier extracted from a fixed character range of \code{path}
#' (\code{str_sub(path, 22, -6)}) — the caller is responsible for ensuring
#' \code{path} has a consistent directory-depth so this offset lines up with
#' the filename portion of the path.
#'
#' @param path Path to a raw LI-6800 \code{.xlsx} export
#' @param type Curve type: one of \code{"A-Ci"}/\code{"ACi"} or
#'   \code{"R-T"}/\code{"RT"}
#' @param area_correction Manually measured leaf area (cm^2) to re-apply via
#'   \code{\link{correct_6800_area}}, or \code{NULL}/\code{NA} to use the
#'   file's own embedded leaf area unmodified
#'
#' @return A tidied data frame with a fixed set of ESS-DIVE-style columns, or
#'   \code{NULL} (with a warning) if \code{path} has no data rows - callers
#'   using \code{bind_rows()}/\code{lapply()}/\code{mapply()} over multiple
#'   paths can pass the result straight through, since \code{bind_rows()}
#'   silently drops \code{NULL} list elements
#' @export
tidy_gasxchg = function(
  path,
  type = c("A-Ci", "R-T"),
  area_correction = NULL
)
{
  # Read in data
  df <-
    if (is.null(area_correction) | is.na(area_correction)) {
      read_6800_xlsx(path)
    } else {
      correct_6800_area(path, area_correction)
    }

  # Some raw files have no data rows (e.g. a curve that was started and
  # immediately aborted); skip them rather than erroring, so a single empty
  # file doesn't take down a whole batch of bind_rows()/mapply() calls
  if (nrow(df) == 0) {
    warning("Skipping ", path, ": file contains no data rows", call. = FALSE)
    return(NULL)
  }

  # Filter out replicated columns, and metadata columns that collide with
  # another column once group prefixes are stripped below (e.g. both
  # TGA_iso.Version and gasanalyzer.Version become "Version"). any_of()
  # rather than a plain vector so a file missing an optional module's
  # columns (e.g. no TGA_iso attachment) doesn't error out here.
  df <- df |>
    select(!any_of(c(
      "SysConst.AvgTime","GasEx.Time", "Meas.Time","MchEvent.Time",
      "MchEvent.HHMMSS",
      "Raw.H2Oa", "Raw.H2Or",
      "TGA_iso.Version", "gasanalyzer.Version"
    )))

  # Remove column categories
  names(df) = gsub(".*\\.","", names(df))

  if (type %in% c("A-Ci", "ACi")) {

    # Get identifying information
    curve_ID = str_sub(path, 22, -6)
    split_ID = str_split(curve_ID, "_", simplify = TRUE)
    sp_code = str_to_upper(split_ID[1])
    elevation = split_ID[2]
    replicate = split_ID[3]
    if (split_ID[4] == "newleaf") {
      replicate = paste(replicate, "newleaf", sep = "_")
    }
    i = length(split_ID)
    meas_temp = split_ID[i]

    df = df |>
      # Add identifying columns
      mutate(
        curve_ID, sp_code, elevation, replicate,
        type = ifelse(split_ID[4] %in% c("25", "25only"), "ACi_25", "ACi-T"),
        meas_temp,
        .after = HHMMSS
      ) |>
      # Rename columns according to ESS-DIVE
      dplyr::rename(
        Patm = Pa,
        RHs = RHcham,
        area = S,
        flow = Flow
      )

  } else if (type %in% c("R-T", "RT")) {

    curve_ID = str_sub(path, 22, -6)
    split_ID = str_split(curve_ID, "_", simplify = TRUE)
    sp_code = str_to_upper(split_ID[1])
    elevation = split_ID[2]
    replicate = split_ID[3]

    df = df |>
      # Add identifying columns
      mutate(
        curve_ID, sp_code, elevation, replicate, type = "RT",
        .after = HHMMSS
      ) |>
      # Rename columns according to ESS-DIVE
      dplyr::rename(
        Patm = Pa,
        RHs = RHcham,
        area = S,
        flow = Flow,
        fan = Fan
      ) |>
      # Correct A to R
      dplyr::mutate(
        Rdark = -A,
        .after = A
      ) |>
      select(!A)
  }

  df = df |>
    select(any_of(cols_to_select))

  # Reformat as separate Date and time columns
  df$Date = ymd_hms(df$Date)
  time = hms::as_hms(df$Date)
  df = df |> add_column(time, .before = 2)

  return(df)
}

#' Add a 4-letter species code
#'
#' Derives a 4-letter uppercase species code (first 2 letters of the genus,
#' first 2 letters of the specific epithet) from a \code{species} column.
#'
#' @param df A data frame with a \code{species} column of \code{"Genus
#'   epithet"} names
#'
#' @return \code{df} with an added \code{sp_code} column, placed after
#'   \code{species}
#' @export
add_sp_code = function(df)
{
  genus = vapply(strsplit(df$species," "), `[`, 1, FUN.VALUE = character(1))
  sp_epithet = vapply(strsplit(df$species," "), `[`, 2, FUN.VALUE = character(1))

  df = df |>
    mutate(sp_code = toupper(paste0(str_sub(genus, 1, 2), str_sub(sp_epithet, 1, 2))),
           .after = species)
  return(df)
}

#' Add species names from a 4-letter species code
#'
#' Maps a \code{sp_code} column back to full species names using a
#' study-specific lookup, rather than a hardcoded table, so it can be reused
#' across studies with different species lists.
#'
#' @param df A data frame with a \code{sp_code} column
#' @param lookup A named character vector mapping \code{sp_code} to species
#'   name (\code{names(lookup)} are the 4-letter codes), or a 2-column data
#'   frame/list with elements \code{sp_code} and \code{species}. Codes not
#'   yet matched to a taxonomic name should still be included with a
#'   placeholder value (e.g. \code{""}) rather than omitted, so they remain
#'   visibly unresolved; codes entirely absent from \code{lookup} map to
#'   \code{NA}.
#'
#' @return \code{df} with an added/overwritten \code{species} column
#' @export
add_species = function(df, lookup)
{
  if (is.data.frame(lookup) || is.list(lookup)) {
    lookup = stats::setNames(lookup[["species"]], lookup[["sp_code"]])
  }

  df = df |>
    mutate(species = unname(lookup[sp_code]))

  return(df)
}
