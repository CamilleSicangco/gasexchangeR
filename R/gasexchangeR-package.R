#' @importFrom dplyr select mutate rename filter group_by summarise arrange relocate bind_rows case_when any_of first across
#' @importFrom rlang sym "!!"
#' @importFrom stringr str_sub
#' @importFrom plantecophys fitacis
#' @importFrom nlme nlme
#' @importFrom lme4 lmer
#' @importFrom stats anova AIC shapiro.test residuals fitted coef lm nls cor as.formula reformulate
#' @importFrom utils globalVariables
#' @importFrom PhotoGEA check_response_curve_data organize_response_curve_data calculate_temperature_response calculate_total_pressure consolidate fit_c3_aci calculate_jmax is.exdf check_required_variables multi_curve_colors

# PhotoGEA's c3_temperature_param_bernacchi and jmax_temperature_param_bernacchi
# are lazy-loaded datasets, and unit_dictionary is an unexported internal
# function - none of the three can be declared via @importFrom, so call sites
# reference them as PhotoGEA::c3_temperature_param_bernacchi /
# PhotoGEA::jmax_temperature_param_bernacchi / PhotoGEA:::unit_dictionary
# instead (see fit_ACis.R).

# These names are all data frame/exdf columns referenced via tidyverse-style
# non-standard evaluation (e.g. inside filter()/mutate()), which R CMD check
# cannot tell apart from genuinely undefined global variables.
utils::globalVariables(c(
  "A", "Ac250", "Aj250", "Ci", "E", "Fan", "Flow", "Gamma_star_tl_avg",
  "HHMMSS", "Jmax_tl", "PARi", "Pa", "Patm", "Photo", "RHcham",
  "RL_tl_avg", "S", "Tirga", "Tleaf", "Tleaf.x", "TleafCnd_avg", "Tp_tl_avg",
  "Vcmax_tl_avg", "curve_ID", "location", "method", "par", "sp_code",
  "species"
))

NULL
