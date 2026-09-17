# Internal PhotoGEA-engine helpers -------------------------------------------

#' Prepare a Licor exdf for A-Ci fitting with PhotoGEA
#'
#' Adds a curve-point-count column, validates/organizes each curve by point
#' count, and computes temperature-dependent C3 photosynthetic parameters and
#' total chamber pressure. Internal helper used by \code{\link{fit_ACis}}
#' when \code{engine = "PhotoGEA"}.
#'
#' @param licor_exdf A Licor \code{exdf} object (see \code{PhotoGEA})
#' @param identifier_columns Character vector of column names identifying an
#'   individual curve
#'
#' @return The prepared \code{exdf}, ready for \code{PhotoGEA::fit_c3_aci()}
#' @noRd
organise_data_for_PhotoGEA = function(
  licor_exdf,
  identifier_columns
)
{
  # Add a new column called `curve_npts` that stores the number of points in each
  # response curve
  licor_exdf <- do.call(rbind, by(licor_exdf, licor_exdf[, identifier_columns], function(x) {
    x[, 'curve_npts'] <- nrow(x)
    x
  }))

  # Check and process each group of curves depending on how many points are in the
  # curve
  gasex_exdf_list_processed <- by(licor_exdf, licor_exdf[, 'curve_npts'], function(x) {
    # Get the number of points in these curves
    npts <- x[1, 'curve_npts']

    # Make sure the data meets basic requirements
    check_response_curve_data(x, identifier_columns, npts, 'Ca', driving_column_tolerance = Inf)

    # Remove points with duplicated `Ca` values and order by `Ca`. Ordering
    # by `Ci` instead (as this used to) is usually equivalent, since Ci
    # normally increases with Ca -- but for any curve with even a local dip
    # or instability in the Ci-vs-Ca relationship, ordering by Ci scrambles
    # the Ca order. PhotoGEA::estimate_operating_point() then interpolates
    # Cc against Ca using `approx(..., ties = list("ordered", mean))`, which
    # requires Ca to already be sorted -- with scrambled Ca, "ordered" tie
    # handling silently returns NA instead of erroring, which cascades into
    # a 0-row exdf assignment and a cryptic "replacement has 1 row, data has
    # 0" error deep inside PhotoGEA. Ordering by Ca directly avoids this;
    # verified this doesn't change fitted Vcmax/Jmax/etc for curves that
    # were already fitting successfully (only the operating-point row
    # selection depends on this order).
    organize_response_curve_data(
      x,
      identifier_columns,
      c(),
      column_for_ordering = 'Ca',
      columns_to_average = c('TleafCnd', 'Qin')
    )
  })

  # Use `rbind` to recombine all the data
  licor_exdf <- do.call(rbind, gasex_exdf_list_processed)

  # Calculate temperature-dependent values of C3 photosynthetic parameters
  # (a lazy-loaded PhotoGEA dataset, so it can't be @importFrom'd). Bernacchi
  # et al. 2001, matching jmax_temperature_param_bernacchi below and every
  # consuming project's own historical PhotoGEA fits -- don't switch to
  # c3_temperature_param_sharkey without checking downstream saved fits, since
  # this changes fitted Vcmax/Km/GammaStar values, not just refactors them.
  licor_exdf <- calculate_temperature_response(licor_exdf, PhotoGEA::c3_temperature_param_bernacchi)

  # Calculate the total pressure in the leaf chamber
  licor_exdf <- calculate_total_pressure(licor_exdf)
  return(licor_exdf)
}

#' Fit A-Ci curves with PhotoGEA
#'
#' Internal helper used by \code{\link{fit_ACis}} when
#' \code{engine = "PhotoGEA"}.
#'
#' @param licor_exdf An \code{exdf} prepared by
#'   \code{\link{organise_data_for_PhotoGEA}}
#' @param identifier_columns Character vector of column names identifying an
#'   individual curve
#' @param Ca_atmospheric Atmospheric CO2 concentration, ppm
#'
#' @return A PhotoGEA fit-result list (\code{fits}, \code{parameters},
#'   \code{fits_interpolated})
#' @noRd
fit_ACis_PhotoGEA = function(
  licor_exdf,
  identifier_columns,
  Ca_atmospheric = 420
)
{
  # PhotoGEA::fit_c3_aci()'s optimizer has a stochastic component; fix the
  # seed so repeated fits (and comparisons against previously saved results)
  # are reproducible.
  set.seed(1234)

  # Fit each curve, catching any per-curve error (e.g. a curve too short/
  # sparse for one of PhotoGEA's internal grid-construction steps) so a
  # single bad curve can't take down the fit for every other curve in the
  # dataset -- by() itself has no such isolation; an uncaught error in any
  # one group aborts the whole call. A failing curve is warned about by
  # name and then dropped from the fit itself (not NA-filled): consolidate()
  # requires every element to share the same list structure, and there's no
  # meaningful partial fit_c3_aci() result to fabricate for a curve whose
  # fit crashed outright. Its identifier is recorded in the "failed_curves"
  # attribute below instead, so callers who want a row per input curve can
  # still recover it via extract_curve_params(..., keep_failed = TRUE).
  failed_curves <- c()
  by_result <- by(
    licor_exdf,                       # The `exdf` object containing the curves
    licor_exdf[, identifier_columns], # A factor used to split `gasex_data` into chunks
    function(x, ...) {
      tryCatch(
        fit_c3_aci(x, ...),
        error = function(e) {
          curve_id <- unique(x[, identifier_columns])
          warning(
            "fit_c3_aci() failed for curve(s) ", paste(curve_id, collapse = ", "),
            ": ", conditionMessage(e), " -- excluded from the fit results.",
            call. = FALSE
          )
          failed_curves <<- c(failed_curves, curve_id)
          NULL
        }
      )
    },
    Ca_atmospheric = Ca_atmospheric   # Additional argument passed to `fit_c3_aci`
  )
  by_result <- by_result[!vapply(by_result, is.null, logical(1))]

  # Fit the A-Ci curves
  ACi_fits <- consolidate(by_result)

  # Recorded as an attribute (not a 4th list element) so this stays
  # transparent to plot_c3_aci_fit_indiv0()'s `length(fit_results) != 3`
  # structural check and to consolidate()-shaped consumers generally.
  attr(ACi_fits, "failed_curves") <- failed_curves

  # Calculate temperature-dependent values of Jmax-related parameters
  # (a lazy-loaded PhotoGEA dataset, so it can't be @importFrom'd)
  ACi_fits$parameters <- calculate_temperature_response(
    ACi_fits$parameters,
    PhotoGEA::jmax_temperature_param_bernacchi,
    'TleafCnd_avg'
  )

  # Calculate Jmax
  ACi_fits$parameters <- calculate_jmax(
    ACi_fits$parameters, 
    alpha_j_at_25 = 0.24, 
    theta_j_at_25 = 0.85
  )

  return(ACi_fits)
}

#' Extract An/Ac/Aj/Ap at a fixed Ci from a PhotoGEA fit result
#'
#' Internal helper used by \code{\link{extract_curve_params}} when
#' \code{engine = "PhotoGEA"}.
#'
#' @param fit_result A PhotoGEA fit-result list from
#'   \code{\link{fit_ACis_PhotoGEA}}
#' @param Ci Target intercellular CO2 concentration, ppm
#' @param val Which interpolated assimilation rate to extract
#' @param identifier_column Name of the per-curve identifier column in
#'   \code{fit_result$fits_interpolated$main_data} (whatever \code{group_cols}/
#'   \code{byname} was originally used for that curve)
#'
#' @return A named numeric vector, one value per curve
#' @noRd
extract_val_fixedCi = function(
  fit_result,
  Ci = 250,
  val = c("An", "Ac", "Aj", "Ap"),
  identifier_column = "curve_ID"
)
{
  fits.l = split(fit_result$fits_interpolated$main_data, fit_result$fits_interpolated$main_data[[identifier_column]])
  val_out = sapply(fits.l, function(df) {
    # A curve whose fit failed entirely (e.g. no reliable Vcmax) has Ci ==
    # NA for every interpolated point, so which.min() finds no minimum and
    # returns integer(0) -- df[[val]][integer(0)] is then length 0 instead
    # of the length-1 result every other curve produces. Left as-is, that
    # breaks sapply()'s simplification to a plain vector for the *entire*
    # batch (falling back to a list), which only surfaces much later as an
    # opaque "comparison of these types is not implemented" when the list
    # is compared with `<`. Degrade explicitly to NA here instead, so one
    # failed curve can't corrupt every other curve's result type.
    if (all(is.na(df$Ci))) {
      return(NA_real_)
    }
    i = which.min(abs(Ci - df$Ci))
    df[[val]][i]
  })
  return(val_out)
}

#' Ci transition point from a plantecophys fit
#' @param fit A single \code{plantecophys::fitaci} result
#' @return The curve's Ci transition point (unnamed)
#' @noRd
get_Ci_transition = function(fit)
{
  unname(fit$Ci_transition)
}

# Public, engine-aware A-Ci fitting API ---------------------------------------

#' Fit A-Ci curves
#'
#' Fits A-Ci curves using either the \pkg{PhotoGEA} or \pkg{plantecophys}
#' package as the fitting engine, replacing the previously separate
#' PhotoGEA- and plantecophys-specific fitting functions with a single
#' entry point.
#'
#' @param data For \code{engine = "PhotoGEA"}, a Licor \code{exdf} object.
#'   For \code{engine = "plantecophys"}, a plain data frame containing (at
#'   least) \code{A}, \code{Ci}, \code{Tleaf}, a PAR column (see
#'   \code{PAR_col}), and the column(s) named in \code{group_cols}; pre-filter
#'   this to whichever subset of curves you want fitted (e.g. one sample or
#'   one location) before calling.
#' @param group_cols Character vector of column name(s) identifying an
#'   individual curve (e.g. \code{"curve_ID"})
#' @param engine Curve-fitting engine: \code{"PhotoGEA"} (default) or
#'   \code{"plantecophys"}
#' @param Ca_atmospheric Atmospheric CO2 concentration, ppm (PhotoGEA engine
#'   only)
#' @param fitTPU Whether to fit a TPU-limited term (plantecophys engine only;
#'   PhotoGEA's model always considers Ac/Aj/Ap limitation)
#' @param PAR_col Name of the incident-PAR column to use as \code{PARi}
#'   (plantecophys engine only)
#' @param Patm_by Optional column name to split \code{data} by before
#'   fitting, using each split's mean \code{Patm} (plantecophys engine only;
#'   e.g. fit separately per \code{"location"} using each location's mean
#'   atmospheric pressure). If \code{NULL} (default), plantecophys's own
#'   default atmospheric pressure handling is used.
#'
#' @return The fitting engine's native fit-result object: a PhotoGEA fit
#'   list (PhotoGEA engine) or a list of \code{plantecophys::fitaci} results
#'   (plantecophys engine), to be passed into \code{\link{extract_curve_params}}
#'   and \code{\link{get_r2}}
#' @export
fit_ACis = function(
  data,
  group_cols,
  engine = c("PhotoGEA", "plantecophys"),
  Ca_atmospheric = 420,
  fitTPU = FALSE,
  PAR_col = "Qin",
  Patm_by = NULL
)
{
  engine = match.arg(engine)

  if (engine == "PhotoGEA") {
    organised = organise_data_for_PhotoGEA(data, identifier_columns = group_cols)
    fits = fit_ACis_PhotoGEA(organised, identifier_columns = group_cols, Ca_atmospheric = Ca_atmospheric)
  } else {
    df = data.frame(data) |>
      dplyr::rename(Photo = A, PARi = !!sym(PAR_col)) |>
      select(any_of(group_cols), Ci, Photo, Tleaf, PARi)

    if (is.null(Patm_by)) {
      fits = fitacis(df, group = group_cols, fitmethod = "bilinear",
                     id = group_cols, Tcorrect = FALSE, fitTPU = fitTPU)
    } else {
      Patm_bygroup = data |>
        select(any_of(Patm_by), Patm) |>
        group_by(dplyr::across(dplyr::all_of(Patm_by))) |>
        summarise(Patm = mean(Patm), .groups = "drop")

      df.l = split(df, data[[Patm_by]])
      fits = lapply(seq_len(nrow(Patm_bygroup)), function(i) {
        fitacis(df.l[[i]], Patm = Patm_bygroup$Patm[i],
               group = group_cols, fitmethod = "bilinear",
               Tcorrect = FALSE, fitTPU = fitTPU)
      })
      fits = unlist(fits, recursive = FALSE)
    }
  }

  return(fits)
}

#' R-squared of measured vs. modeled assimilation
#'
#' @param fit_result A fit-result object from \code{\link{fit_ACis}}
#' @param engine Curve-fitting engine the fit came from: \code{"PhotoGEA"}
#'   (default) or \code{"plantecophys"}
#' @param identifier_column Name of the per-curve identifier column in
#'   \code{fit_result$fits$main_data} (PhotoGEA engine only; whatever
#'   \code{group_cols} was originally used for that curve)
#'
#' @return A named numeric vector of R^2 values, one per curve
#' @export
get_r2 = function(
  fit_result,
  engine = c("PhotoGEA", "plantecophys"),
  identifier_column = "curve_ID"
)
{
  engine = match.arg(engine)

  if (engine == "PhotoGEA") {
    fits.l = split(fit_result$fits$main_data, fit_result$fits$main_data[[identifier_column]])
    R2 = sapply(fits.l, function(df) {
      lmfit = lm(A ~ A_fit, df)
      summary(lmfit)$r.squared
    })
  } else {
    R2 = sapply(fit_result, function(x) {
      lmfit = lm(Ameas ~ Amodel, data = x$df)
      summary(lmfit)$r.squared
    })
  }

  return(R2)
}

#' Extract fitted A-Ci curve parameters
#'
#' Builds a data frame of fitted Vcmax/Jmax-related quantities (including A
#' at fixed Ci = 250/275/360 ppm, Rubisco- vs. RuBP-regeneration-limitation
#' flags, and R^2), merged back onto the original curve-level data.
#' Replaces the previously separate PhotoGEA- and plantecophys-specific
#' extraction functions.
#'
#' @param fit_result A fit-result object from \code{\link{fit_ACis}}
#' @param data The original curve-level data frame to merge fitted
#'   parameters onto
#' @param byname Column name(s) to merge \code{fit_result}'s parameters onto
#'   \code{data} by. For the PhotoGEA engine this is also the per-curve
#'   identifier column used internally for R^2 and fixed-Ci extraction (the
#'   same column \code{group_cols} named when calling \code{\link{fit_ACis}}),
#'   not a separate assumption of its own.
#' @param engine Curve-fitting engine the fit came from: \code{"PhotoGEA"}
#'   (default) or \code{"plantecophys"}
#' @param Ci_targets Intercellular CO2 concentrations (ppm) to extract
#'   A/Ac/Aj/Ap at (PhotoGEA engine only; ignored for \code{"plantecophys"}).
#'   Defaults to \code{c(250, 275, 360)}, matching this function's original,
#'   fixed set of extraction points.
#' @param average_cols Length-2 character vector naming the two ends of a
#'   positional range of columns in the merged data (PhotoGEA engine only;
#'   ignored for \code{"plantecophys"}) to average across a curve's points
#'   when collapsing to one row per curve, rather than just keep the first
#'   value like every other merged column -- e.g. raw per-point gas-exchange
#'   readings that vary within a curve. Defaults to \code{c("E", "Tirga")};
#'   pass the two columns bounding the equivalent range in your own data if
#'   they're named differently (e.g. \code{c("gsw", "Pa")}).
#' @param keep_failed Whether to keep curves \code{fit_ACis()} failed to fit
#'   (PhotoGEA engine only; ignored for \code{"plantecophys"}) as extra rows
#'   with all fitted-parameter columns \code{NA} and \code{fit_failed = TRUE}
#'   (successfully-fit rows get \code{fit_failed = FALSE}), instead of the
#'   default of omitting them entirely. Curve-level columns from \code{data}
#'   itself (e.g. \code{sp_code}) are still populated on these rows.
#'
#' @return \code{data}, with fitted curve parameters merged on
#' @export
extract_curve_params = function(
  fit_result,
  data,
  byname,
  engine = c("PhotoGEA", "plantecophys"),
  Ci_targets = c(250, 275, 360),
  average_cols = c("E", "Tirga"),
  keep_failed = FALSE
)
{
  engine = match.arg(engine)

  if (engine == "PhotoGEA") {
    extract_curve_params_photogea(fit_result, data, byname, Ci_targets = Ci_targets,
                                   average_cols = average_cols, keep_failed = keep_failed)
  } else {
    extract_curve_params_plantecophys(fit_result, data, byname)
  }
}

#' @noRd
extract_curve_params_plantecophys = function(fit, data, byname)
{
  ret <- bind_rows(sapply(fit, coef, simplify = FALSE), .id = "curve_ID")
  ret$Citrans <- sapply(fit, get_Ci_transition)
  ret$Jmax[ret$Jmax > 1000] <- NA
  ret$rmse <- sapply(fit, "[[", "RMSE")
  ret$R2 <- get_r2(fit, engine = "plantecophys")

  # extract A at fixed Ci, and average Tleaf
  ret$Ac250 <- sapply(fit, function(x) x$Photosyn(Ci = 250)$Ac)
  ret$Aj250 <- sapply(fit, function(x) x$Photosyn(Ci = 250)$Aj)
  ret$Jlim250 <- (ret$Aj250 < ret$Ac250)
  ret$A250 <- sapply(fit, function(x) min(x$Photosyn(Ci = 250)$Ac, x$Photosyn(Ci = 250)$Aj))
  ret$Ac275 <- sapply(fit, function(x) x$Photosyn(Ci = 275)$Ac)
  ret$Aj275 <- sapply(fit, function(x) x$Photosyn(Ci = 275)$Aj)
  ret$Jlim275 <- (ret$Aj275 < ret$Ac275)
  ret$A275_old <- sapply(fit, function(x) min(x$Photosyn(Ci = 275)$Ac, x$Photosyn(Ci = 275)$Aj))
  ret$A275 <- sapply(fit, function(x) x$Photosyn(Ci = 275)$ALEAF)
  ret$Ac360 <- sapply(fit, function(x) x$Photosyn(Ci = 360)$Ac)
  ret$Aj360 <- sapply(fit, function(x) x$Photosyn(Ci = 360)$Aj)
  ret$Jlim360 <- (ret$Aj360 < ret$Ac360)
  ret$A360 <- sapply(fit, function(x) min(x$Photosyn(Ci = 360)$Ac, x$Photosyn(Ci = 360)$Aj))
  ret$Tleaf <- sapply(fit, function(x) mean(x$df$Tleaf))
  ret$TsK <- ret$Tleaf + 273.15
  ret$GammaStar <- sapply(fit, function(x) x$GammaStar)
  ret$Km <- sapply(fit, function(x) x$Km)

  # merge with original data frame to get info on curve
  ret <- merge(ret, data, by = byname)
  ret <- ret |> dplyr::relocate(location:method, .after = curve_ID)
  ret <- ret[!duplicated(ret[, c(byname)]), ]
  ret <- ret |> dplyr::rename(Tleaf = Tleaf.x)

  return(ret)
}

#' @noRd
extract_curve_params_photogea = function(fits_exdf, data, byname, Ci_targets = c(250, 275, 360),
                                          average_cols = c("E", "Tirga"), keep_failed = FALSE)
{
  ret = fits_exdf$parameters$main_data

  # Determine the merge keys now, from the untouched pass-through columns,
  # before any renaming or derived columns are added below. Otherwise a
  # later-created column name that happens to collide with an unrelated
  # column already in `data` (e.g. renaming TleafCnd_avg to Tleaf below,
  # when `data` already has its own, differently-grained, per-point Tleaf
  # column) gets mistaken for a shared identifier and pulled into the
  # merge, silently dropping every row since a single per-curve average
  # essentially never equals any individual point's raw value.
  col_match = intersect(colnames(data), colnames(ret))

  ret$Jmax_tl[is.infinite(ret$Jmax_tl)] <- NA
  ret$R2 <- get_r2(fits_exdf, engine = "PhotoGEA", identifier_column = byname)

  # extract A/Ac/Aj/Ap at each fixed Ci target, and flag Rubisco- vs.
  # RuBP-regeneration-limitation, for every value in Ci_targets (was three
  # hardcoded 250/275/360 blocks; generalized so callers needing an
  # additional extraction point, e.g. 300, don't need a package change)
  for (ci in Ci_targets) {
    ret[[paste0("Ac", ci)]] = extract_val_fixedCi(fits_exdf, Ci = ci, val = "Ac", identifier_column = byname)
    ret[[paste0("Aj", ci)]] = extract_val_fixedCi(fits_exdf, Ci = ci, val = "Aj", identifier_column = byname)
    ret[[paste0("Ap", ci)]] = extract_val_fixedCi(fits_exdf, Ci = ci, val = "Ap", identifier_column = byname)
    ret[[paste0("Jlim", ci)]] = (ret[[paste0("Aj", ci)]] < ret[[paste0("Ac", ci)]] &
                                    ret[[paste0("Aj", ci)]] < ret[[paste0("Ap", ci)]])
    ret[[paste0("A", ci)]] = extract_val_fixedCi(fits_exdf, Ci = ci, val = "An", identifier_column = byname)
  }

  ret = ret |>
    dplyr::rename(Tleaf = TleafCnd_avg,
                  Vcmax = Vcmax_tl_avg,
                  Jmax = Jmax_tl,
                  Rd = RL_tl_avg,
                  TPU = Tp_tl_avg,
                  GammaStar = Gamma_star_tl_avg) |>
    mutate(TsK = Tleaf + 273.15)

  # merge with original data frame to get info on curve
  ret = merge(ret, data, by = col_match)

  # If `data` has its own Tleaf column (a per-point raw reading, distinct
  # from the fitted curve-level average just renamed above), it isn't a
  # merge key -- by design, see col_match above -- so merge() suffixes both
  # as Tleaf.x/Tleaf.y instead of colliding. Keep the fitted average as the
  # curve-level `Tleaf` (matching extract_curve_params_plantecophys()'s
  # equivalent Tleaf.x/Tleaf.y handling) and drop the raw duplicate.
  if ("Tleaf.x" %in% colnames(ret)) {
    ret = ret |>
      dplyr::rename(Tleaf = Tleaf.x) |>
      dplyr::select(-Tleaf.y)
  }

  # average_cols names the two ends of a positional range of per-point
  # columns (e.g. raw gas-exchange readings) to average across a curve's
  # points, rather than keep just the first value like every other column;
  # resolved against ret's actual column order at this point rather than
  # hardcoded, since which columns fall "inside the range" depends on where
  # these two boundary columns happen to sit in a given caller's own data.
  range_idx = which(names(ret) == average_cols[1]):which(names(ret) == average_cols[2])
  range_cols = names(ret)[range_idx]

  ret = ret |>
    group_by(!!sym(byname)) |>
    summarise(across(!dplyr::all_of(range_cols), first),
              across(dplyr::all_of(range_cols), mean),
              .groups = "drop")

  # Re-attach curves fit_ACis() couldn't fit as extra rows, if asked to. Their
  # fitted-parameter columns (Vcmax, R2, A250, ...) are absent from `data`, so
  # bind_rows() naturally leaves them NA -- only the curve-level columns
  # already present in `data` itself (e.g. sp_code) get filled in.
  if (keep_failed) {
    failed = attr(fits_exdf, "failed_curves")

    ret$fit_failed = FALSE

    if (length(failed) > 0) {
      failed_rows = data[data[[byname]] %in% failed, ] |>
        group_by(!!sym(byname)) |>
        summarise(across(everything(), first), .groups = "drop") |>
        mutate(fit_failed = TRUE)

      ret = dplyr::bind_rows(ret, failed_rows)
    }
  }

  return(ret)
}

# A-Ci fit diagnostic plots (PhotoGEA only) -----------------------------------

#' Plot a single PhotoGEA A-Ci fit
#'
#' Diagnostic plot of one curve's measured points and fitted An/Ac/Aj/Ap/Ad
#' lines, with the estimated operating point marked.
#'
#' @param fit_results A PhotoGEA fit-result list from
#'   \code{\link{fit_ACis}}
#' @param identifier_column_name Name of the curve identifier column
#' @param identifier_column_val Curve identifier value to plot
#'
#' @return A \pkg{lattice} plot object
#' @export
plot_c3_aci_fit_indiv = function(
  fit_results,
  identifier_column_name,
  identifier_column_val
)
{
  # Choose line settings
  assim_cols <- multi_curve_colors()[1:5]
  assim_cols[1] <- '#676767'

  line_settings <- list(
    col = assim_cols,
    lwd = c(4, 2, 2, 2, 2),
    lty = c(1, 5, 5, 5, 2)
  )

  # Look this curve's row up by value, not position -- identifier_column_val
  # is the curve's identifier *value* (e.g. a curve_ID string), not a row
  # number, so indexing main_data's columns directly with it (as this used
  # to) did name-based indexing into an unnamed vector and always resolved
  # to NA/blank. Filtering the exdf first and taking row 1 matches how
  # curve_interpolated/curve_data are already selected below.
  curve_parameters <-
    fit_results$parameters[fit_results$parameters[, identifier_column_name] == identifier_column_val, ]

  title = paste0(curve_parameters$main_data$curve_ID[1], ": ",
                curve_parameters$main_data$location[1], ", ",
                curve_parameters$main_data$sample[1], ", Tleaf = ",
                round(curve_parameters$main_data$TleafCnd_avg[1], digits = 1))

  curve_interpolated <-
    subset(fit_results$fits_interpolated$main_data, curve_ID == identifier_column_val)

  # Scale the y-axis to the fitted An curve rather than letting lattice's
  # default auto-scaling span Ac/Aj/Ap/Ad too -- Ap in particular can shoot
  # off to a near-vertical asymptote at low Ci, which would otherwise
  # squash the An/Ac/Aj curves into an unreadable sliver.
  ymin <- min(curve_interpolated$An) - 2
  ymax <- max(curve_interpolated$An) + 5

  lattice::xyplot(An + Ac + Aj + Ap + Ad ~ Ci,
                  curve_interpolated,
                  type = 'l',
                  par.settings = list(superpose.line = line_settings),
                  auto.key = list(space = 'right', lines = TRUE, points = FALSE),
                  ylab = expression(A[n]~(mu*mol~m^-2~s^-1)),
                  xlab = expression(C[i]*" (ppm)"),
                  main = title,
                  ylim = c(ymin, ymax),
                  panel = function(...) {
                    curve_data <-
                      fit_results$fits[fit_results$fits[, identifier_column_name] == identifier_column_val, ]

                    used_for_fit <- points_for_fitting(curve_data)

                    # Plot the fit lines
                    lattice::panel.xyplot(...)

                    # Plot the measured data points
                    lattice::panel.points(
                      curve_data[used_for_fit, "A"] ~ curve_data[used_for_fit, "Ci"],
                      col = 'black',
                      pch = 16
                    )

                    lattice::panel.points(
                      curve_data[!used_for_fit, "A"] ~ curve_data[!used_for_fit, "Ci"],
                      col = 'black',
                      pch = 1
                    )

                    # Plot the operating point
                    lattice::panel.points(
                      curve_parameters[1, 'operating_An_model'] ~ curve_parameters[1, 'operating_Ci'],
                      col = 'red',
                      pch = 1
                    )
                  }
  )
}

#' Plot PhotoGEA A-Ci fits, faceted by curve
#'
#' Like \code{\link{plot_c3_aci_fit_indiv}}, but plots every curve in
#' \code{fit_results} at once, faceted by \code{identifier_column_name}, and
#' supports plotting against either Ci or Cc.
#'
#' @param fit_results A PhotoGEA fit-result list (must have \code{fits},
#'   \code{parameters}, and \code{fits_interpolated} \code{exdf} elements)
#' @param identifier_column_name Name of the curve identifier column
#' @param x_name Name of the x-axis column: must match \code{ci_column_name}
#'   or \code{cc_column_name}
#' @param plot_operating_point Whether to mark the estimated operating point
#' @param plot_Ad Whether to also plot the `Ad` (dark respiration) line
#' @param a_column_name Name of the measured assimilation column
#' @param cc_column_name Name of the chloroplast CO2 column
#' @param ci_column_name Name of the intercellular CO2 column
#' @param ... Additional arguments passed to \code{lattice::xyplot}
#'
#' @return A \pkg{lattice} plot object
#' @export
plot_c3_aci_fit_indiv0 <- function(
  fit_results,
  identifier_column_name,
  x_name,
  plot_operating_point = TRUE,
  plot_Ad = FALSE,
  a_column_name = 'A',
  cc_column_name = 'Cc',
  ci_column_name = 'Ci',
  ...
)
{
  if (length(fit_results) != 3 || !is.exdf(fit_results$fits) ||
      !is.exdf(fit_results$parameters) || !is.exdf(fit_results$fits_interpolated)) {
    stop('fit_results must be the output from fit_c3_aci or fit_c3_variable_j')
  }

  # Make sure the x variable is acceptable
  if (!x_name %in% c(ci_column_name, cc_column_name)) {
    stop('The x_name must be the Ci or Cc column name')
  }

  # Get the appropriate operating point name
  operating_x <- if (x_name == ci_column_name) {
    'operating_Ci'
  } else {
    'operating_Cc'
  }

  # Make sure the required variables are defined and have the correct units
  required_variables <- list()
  required_variables[[cc_column_name]]         <- 'micromol mol^(-1)'
  required_variables[[ci_column_name]]         <- 'micromol mol^(-1)'
  required_variables[[identifier_column_name]] <- NA

  # Don't throw an error if some columns are all NA
  check_required_variables(fit_results$fits_interpolated, required_variables, check_NA = FALSE)

  # unit_dictionary() is an unexported PhotoGEA internal, hence `:::`
  required_variables[[a_column_name]] <- PhotoGEA:::unit_dictionary('A')

  # Don't throw an error if some columns are all NA
  check_required_variables(fit_results$fits, required_variables, check_NA = FALSE)

  required_variables <- list()
  required_variables[['operating_Ci']]       <- 'micromol mol^(-1)'
  required_variables[['operating_Cc']]       <- 'micromol mol^(-1)'
  required_variables[['operating_An_model']] <- PhotoGEA:::unit_dictionary('A')

  # Don't throw an error if some columns are all NA
  check_required_variables(fit_results$parameters, required_variables, check_NA = FALSE)

  # Choose line settings
  assim_cols <- multi_curve_colors()[1:5]
  assim_cols[1] <- '#676767'

  line_settings <- list(
    col = assim_cols,
    lwd = c(4, 2, 2, 2, 2),
    lty = c(1, 5, 5, 5, 2)
  )

  # Specify the y-axis label
  y_label <- paste0(
    'Net CO2 assimilation rate [ ', fit_results$fits_interpolated$units[['An']],
    ' ]\n(filled black circles: measured data used for fits',
    '\nopen black circles: measured data excluded from fits'
  )

  y_label <- paste0(
    y_label,
    if (plot_operating_point) {
      '\nopen red circle: estimated operating point)'
    } else {
      ')'
    }
  )

  # Plot the fits, operating point, and raw data
  lattice::xyplot(
    if (plot_Ad) {
      An + Ac + Aj + Ap + Ad ~ fit_results$fits_interpolated[, x_name] | fit_results$fits_interpolated[, identifier_column_name]
    } else {
      An + Ac + Aj + Ap ~ fit_results$fits_interpolated[, x_name] | fit_results$fits_interpolated[, identifier_column_name]
    },
    data = fit_results$fits_interpolated$main_data,
    type = 'l',
    par.settings = list(superpose.line = line_settings),
    auto.key = list(space = 'right', lines = TRUE, points = FALSE),
    xlab = paste(x_name, '[', fit_results$fits_interpolated$units[[x_name]], ']'),
    ylab = y_label,
    curve_ids = fit_results$fits_interpolated[, identifier_column_name],
    panel = function(...) {
      # Get info about this curve
      args <- list(...)
      curve_id <- args$curve_ids[args$subscripts][1]

      curve_parameters <-
        fit_results$parameters[fit_results$parameters[, identifier_column_name] == curve_id, ]

      curve_data <-
        fit_results$fits[fit_results$fits[, identifier_column_name] == curve_id, ]

      used_for_fit <- points_for_fitting(curve_data)

      # Plot the fit lines
      lattice::panel.xyplot(...)

      # Plot the measured data points
      lattice::panel.points(
        curve_data[used_for_fit, a_column_name] ~ curve_data[used_for_fit, x_name],
        col = 'black',
        pch = 16
      )

      lattice::panel.points(
        curve_data[!used_for_fit, a_column_name] ~ curve_data[!used_for_fit, x_name],
        col = 'black',
        pch = 1
      )

      # Plot the operating point, if desired
      if (plot_operating_point) {
        lattice::panel.points(
          curve_parameters[1, 'operating_An_model'] ~ curve_parameters[1, operating_x],
          col = 'red',
          pch = 1
        )
      }
    },
    ...
  )
}

#' Identify points included when fitting a curve
#'
#' If \code{curve_data} has an \code{include_when_fitting} logical column,
#' that's used to determine which points were included; otherwise all points
#' are assumed included.
#'
#' @param curve_data A curve's fit-result data frame/exdf
#'
#' @return A logical vector, one per row of \code{curve_data}
#' @noRd
points_for_fitting <- function(curve_data)
{
  if ('include_when_fitting' %in% colnames(curve_data)) {
    as.logical(curve_data[, 'include_when_fitting'])
  } else {
    rep_len(TRUE, nrow(curve_data))
  }
}
