#' Response (y) Statistics Within Levels/Bins of a Predictor (x)
#' 
#' Statistics for \code{y} are plotted with respect to each level or bin of \code{x}. Plotted 
#' statistics can be proportions, log-odds, or weight-of-evidence values. Bins can be created 
#' using raw factor levels, quantile breakpoints, uniform breakpoints, or recursive partitioning. 
#' Additional arguments may be passed to \code{rpart.control()} to fine-tune recursive partitioning. 
#' Plots showing the \code{ymetric} for each value of \code{xsplit} as well as the total volume in 
#' each bin are printed to the current graphics device. In addition, two measures of the overall 
#' strength of the predictive relationship (Information Value & ChiSq) are calculated and returned.
#'
#' @param y (numeric) binary response vector
#' @param x (numeric) numeric or factor predictor vector
#' @param ymetric (character) statistic to calculate for \code{y}: \code{c('proportion', 'logodds', 'woe')}
#' @param xsplit (character) method used to bin \code{x}: \code{c('quantile', 'uniform', 'rpart')}
#' @param nbins (numeric) number of bins to create from \code{x}
#' @param nabin (logical) whether to include an additional bin for missing \code{x} values
#' @param yticks (numeric) number of tick marks to display on the y-axis of plots
#' @param ... (args) additional arguments to pass to \code{rpart.control()}
#' 
#' @details 
#' If \code{xsplit='rpart'} bins will be created based on recursive partitioning for both 
#' numeric and factor variables and the \code{nbins} argument will be ignored. Pass additional 
#' control parameters (e.g. cp, minbucket) in the function call to control partitioning 
#' behavior. If zero or greater than 20 bins are created using the rpart control settings 
#' passed the function will throw an error. If x is a factor variable the x-axis labels on 
#' the returned plots will correspond to the index positions of the levels of x (and not 
#' the factor labels themselves) in each bin. It's generally not a good idea to use recursive 
#' partitioning with more than 50 factor levels. If x is a numeric variable the x-axis labels 
#' will be the range cutpoints for each bin created via recursive partitioning.
#' 
#' If \code{xsplit=c('uniform','quantile')} and x is a factor variable its levels are used 
#' directly as bins and the \code{nbins} argument will be ignored. If x is a numeric variable 
#' bins are calculated by dividing the range of x into buckets of either equal size (uniform) 
#' or equal count (quantile). If quantile breakpoints are not unique then adjacent identical 
#' bins will be combined. 
#' 
#' If bins get created which have either zero volume or zero variance then log-odds and woe
#' cannot be calculated. Any such bins will be excluded from both the displayed plots and 
#' also the calculation of information value for the variable. This problem can typically 
#' be solved by using quantile binning and/or reducing the number of bins created.
#' 
#' @return a list containing the following elements:
#' \item{iv} {Information Value}
#' \item{chi2} {ChiSq Statistic} 
#' \item{yPlot} {ggplot object of \code{ymetric} vs. \code{bins}}
#' \item{vPlot} {ggplot object of bin sizes or volume}
#'  
#' @examples 
#' data(diamonds, package = 'ggplot2')
#' y  <- as.numeric(diamonds$price > mean(diamonds$price))
#' x1 <- diamonds$carat
#' x2 <- diamonds$clarity
#' x3 <- diamonds$y
#' 
#' res <- plotYXbin(y, x1) 
#' res <- plotYXbin(y, x1, nbins = 8, nabin = FALSE)
#' res <- plotYXbin(y, x2, ymetric = 'woe')
#' res <- plotYXbin(y, x3, ymetric = 'proportion', xsplit = 'rpart', cp = 1e-4, minbucket = 100)
#' @export

plotYXbin <- function(y, x, ymetric='proportion', xsplit='quantile', nbins=10, nabin=TRUE, yticks=6, ...) {
    
    # require external packages
    
    if (!require(ggplot2) || !require(rpart) || !require(grid) || !require(gridExtra)) {
        stop('packages [ggplot2, magrittr, rpart, grid, gridExtra] are required for this function')
    }
    
    # parameter error handling (y, x)
    
    if ((length(unique(y)) != 2) || (min(y) != 0) || (max(y) != 1)) {
        stop('y must be a binary numeric variable with no missing values')
    }
    if (!(typeof(x) %in% c('integer', 'double'))) {
        stop('x must be a numeric or factor variable')
    }
    
    # parameter handling (xsplit) and create x bins
    
    if (xsplit == 'rpart') {
        tree.fit <- rpart::rpart(y ~ x, method = 'class', control = rpart::rpart.control(...)) 
        if (length(tree.fit$splits) == 0) {
            stop('no splits produced with current rpart control settings')
        }
        if (length(tree.fit$frame$var[tree.fit$frame$var == '<leaf>']) > 20) {
            stop('too many leaf nodes (>20) produced with current rpart control settings')
        }
        if (is.factor(x)) {
            letvec  <- c(letters, LETTERS)
            isplits <- length(tree.fit$frame$var[tree.fit$frame$var == 'x'])
            clabels <- strsplit(gsub('^x=', '', labels(tree.fit)[(isplits+1):length(labels(tree.fit))]), '')
            ilabels <- sapply(clabels, function(x) paste(sapply(strsplit(x, ''), function(x) paste(which(letvec==x))), collapse=','))
            bins    <- factor(tree.fit$where, labels = ilabels)
            nbins   <- length(levels(bins))
        }
        else {
            splits <- c(min(x, na.rm=T), sort(tree.fit$splits[,'index']), max(x, na.rm=T))
            bins   <- cut(x, splits, right = FALSE, include.lowest = TRUE)
            nbins  <- length(levels(bins))
        }
    }
    else if (is.factor(x) || length(unique(x)) <= nbins) {
        bins  <- as.factor(x)
        nbins <- length(levels(bins))
    }
    else if (xsplit == 'uniform') {
        bins <- cut(x, nbins, include.lowest = TRUE)
    }
    else if (xsplit == 'quantile') {
        bins  <- cut(x, unique(quantile(x, seq(0, 1, 1/nbins), na.rm = TRUE)), include.lowest = TRUE)
        nbins <- length(levels(bins))
    }
    else {
        stop('xsplit must be one of [quantile, uniform, rpart]')
    }
    
    # include a bin for missing values if any exist
    
    if (nabin == TRUE && length(x[is.na(x)])>0) {
        levels(bins)[nbins+1] <- 'NA'
        bins[is.na(bins)]     <- 'NA'
        nbins <- nbins + 1
    }
    
    # parameter processing (ymetric)
    
    if (ymetric == 'proportion') {
        yfun <- function(y) mean(y, na.rm=T)
        ylab <- 'Proportion (Y)'
    }
    else if (ymetric == 'logodds') {
        yfun <- function(y) ifelse((mean(y, na.rm=T)>0 && mean(y, na.rm=T)<1), log(mean(y, na.rm=T)/(1-mean(y, na.rm=T))), NA)
        ylab <- 'Log-Odds (Y)'
    }
    else if (ymetric == 'woe') {
        ylab <- 'WOE (X)'
    }
    else {
        stop('ymetric must be one of [proportion, logodds, woe]')
    }
    
    # calculate Information Value (IV) and ChiSq statistics
    
    px_y0    <- table(bins[y==0])/length(bins[y==0])
    px_y1    <- table(bins[y==1])/length(bins[y==1])
    woe      <- ifelse(px_y0 > 0 & px_y1 > 0, log(px_y1/px_y0), NA)
    iv       <- sum((px_y1 - px_y0) * woe, na.rm = TRUE)
    chi2     <- suppressWarnings(chisq.test(table(bins, y)))$statistic
    
    # calculate (x, y) values to plot
    
    xvalues <- 1:nbins
    xlabels <- levels(bins)
    xfreqs  <- table(bins)
    
    if (ymetric == 'woe') {
        yvalues <- woe
    }
    else {
        yvalues <- tapply(y, bins, yfun)
    }
    
    ylimits <- c(min(yvalues, na.rm=T), max(yvalues, na.rm=T))
    ybuffer <- 0.1 * (ylimits[2] - ylimits[1])
    ybreaks <- seq(ylimits[1] - ybuffer, ylimits[2] + ybuffer, length.out = yticks)
    vbreaks <- seq(0, max(xfreqs), length.out = yticks)
    
    # produce the (x, y) and x-volume plots
    
    yPlot <- ggplot2::ggplot(data.frame(yval = yvalues, xval = xvalues)) + 
    geom_bar(aes(x = xval, y = yval), width=0.95, alpha=0.75, color='black', fill='blue', stat='identity') +
    geom_line(aes(x = xval, y = yval), size = 1, color = '#ff3399') +
    geom_point(aes(x = xval, y = yval), fill = 'black', size = 2) +
    coord_cartesian(ylim = c(min(ybreaks), max(ybreaks))) +
    scale_x_continuous(name = 'Bin Values/Ranges (X)', breaks = 1:nbins, labels = xlabels) +
    scale_y_continuous(name = ylab, breaks = round(ybreaks, 2))
    
    vPlot <- ggplot2::ggplot(data.frame(xval = as.numeric(bins[!is.na(bins)]))) + 
    geom_bar(aes(x = xval), width=0.95, alpha=0.75, stat='count', color='black', fill='purple') +
    scale_x_continuous(name = 'Bin Values/Ranges (X)', breaks = 1:nbins, labels = xlabels) +
    scale_y_continuous(name = 'Volume (X)', breaks = vbreaks, labels = as.character(round(vbreaks, -1))) +
    labs(x = 'Bin Values/Ranges (X)')
    
    # print both plots and return all results invisibly in a list
    
    suppressWarnings(grid.arrange(yPlot, vPlot))
    print(paste("Information Value:", round(iv, 3), "|", "ChiSq Statistic:", round(chi2, 2)))
    invisible(list(iv = iv, chi2 = chi2, yPlot = yPlot, vPlot = vPlot))
    
}
