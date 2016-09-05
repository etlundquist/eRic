#' Calculate summary statistics for time series variables
#' For each combination of [id/measure] summary statistics are calculated over all periods
#'
#' @param df (data frame) data set containing the time series variables 
#' @param idvar (character) name of the variable that uniquely identifies subjects in the data
#' @param pfx (character) vector of prefix patterns for each time series variable
#' @param sfx (character) pattern yielding period-specific values for each time series variable
#' @param ... (functions) functions to evaluate for each [id/measure] (built-in or user-defined)
#' @return a data frame containing the \code{idvar} and one column for each combination of [measure/stat]
#' 
#' @importFrom magrittr %>%
#' @export

tsSummary <- function(df, idvar, pfx, sfx, ...) {
    
    # set up regular expressions to extract prefix and time period information
    pregex <- paste0('(', paste(pfx, collapse = '|'), ')')
    cregex <- paste0(pregex, sfx)
    
    # group the data by [id/measure] and evaluate all passed functions on groups
    dplyr::select(df, matches(idvar), matches(cregex)) %>%
    tidyr::gather_('variable', 'value', names(df)[stringr::str_detect(names(df), cregex)]) %>%
    dplyr::mutate(measure = stringr::str_extract(variable, pregex)) %>%
    dplyr::group_by_(idvar, 'measure') %>% 
    dplyr::summarize_each(funs(...), value) %>%
    tidyr::gather(stat, value, -(id:measure)) %>%
    dplyr::mutate(field = paste(measure, stat, sep = '_'), measure = NULL, stat = NULL) %>%
    tidyr::spread(field, value)
    
}
