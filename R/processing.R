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
