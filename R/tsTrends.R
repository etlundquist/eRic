#' Calculate linear trends for a set of time series variables
#' For each combination of \code{idvar} and \code{measure} regress measure values against time period.
#' The result will be the best linear approximation of the measure slope as a function of time
#'
#' @param df (data frame) data frame containing the time series variables of interest
#' @param idvar (character) name of the variable that uniquely identifies observations in the data
#' @param pfx (character) vector of prefix patterns for each time series measure
#' @param sfx (character) pattern to match time period suffix values
#' @param reverse (logical) reverse the order of period values when calculating linear trends
#' @return a data frame containing the \code{idvar} and one column for each time series measure containing the 
#' linear trends for that \code{idvar} and \code{measure} combination
#' 
#' @importFrom magrittr %>%
#' @export

tsTrends <- function(df, idvar, pfx, sfx, reverse = FALSE) {
    
    # set up regular expressions to extract measure prefixes and time period suffixes
    pregex <- paste0('(', paste(pfx, collapse = '|'), ')')
    sregex <- sfx
    cregex <- paste0(pregex, sregex)
    
    # reshape the data long to the [id/measure/time] level
    ldata <- dplyr::select(df, matches(idvar), matches(cregex)) %>%
             tidyr::gather_('variable', 'value', names(df)[stringr::str_detect(names(df), cregex)]) %>%
             dplyr::mutate(measure = stringr::str_extract(variable, pregex), 
                           time    = stringr::str_extract(variable, sregex))
    
    # create an integer time period variable and reverse original suffix sort order if specified
    if (reverse == TRUE) {
        ldata <- dplyr::group_by_(ldata, idvar, 'measure') %>% 
                 dplyr::arrange(desc(time)) %>% dplyr::mutate(time = seq_along(time))
    }
    else {
        ldata <- dplyr::group_by_(ldata, idvar, 'measure') %>% 
                 dplyr::arrange(time) %>% dplyr::mutate(time = seq_along(time))
    }
  
    # regress measure value on time for each [id/measure] and return the slope coefficient
    dplyr::do(ldata, data.frame(slope = lm(value ~ time, data = .)$coefficients[['time']])) %>% 
    dplyr::ungroup() %>% tidyr::spread_('measure', 'slope')
    
}
                    