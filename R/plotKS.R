#' KS Plot & Statistic
#' 
#' Plot the cumulative distributions of each class by a score variable and calculate the KS
#' statistic denoting the max vertical distance between the two cumulative distributions.
#'
#' @param a (numeric) vector of scores for the first class 
#' @param b (numeric) vector of scores for the second class
#' @param la (character) legend label for the first class 
#' @param lb (character) legend label for the second class 
#' 
#' @return a list containing the following elements:
#' \item{ks.stat} {KS statistic}
#' \item{ks.plot} {a ggplot2 object for the KS plot} 
#'  
#' @examples 
#' a <- rnorm(5000, 0, 2)
#' b <- rnorm(5000, 1, 3)
#' res <- plotKS(a, b, la = 'Non-Defaults', lb = 'Defaults')
#' res$ks.stat
#' res$ks.plot
#' ks.test(a, b)
#' @export

plotKS <- function(a, b, la = 'Distribution A', lb = 'Distribution B') {
    
    if (!require(ggplot2)) {
        stop('this function requires package: ggplot2')
    }
    if (!is.vector(a, mode = 'numeric') || !is.vector(b, mode = 'numeric')) {
        stop('a and b must be numeric vectors')
    }
    if (!is.character(la) || !is.character(lb)) {
        stop('la and lb must be character labels')
    }
    if (length(a) < 100 || length(b) < 100) {
        warning('CDF comparisons are unstable using distributions with so few observed values')
    }
    
    cdf.a <- ecdf(a)
    cdf.b <- ecdf(b)
    
    x.range <- seq(min(c(a, b)), max(c(a, b)), length.out = length(c(a, b)))
    ks.stat <- max(abs(cdf.a(x.range) - cdf.b(x.range)))
    
    x.star <- x.range[abs(cdf.a(x.range) - cdf.b(x.range)) == ks.stat][1]
    a.star <- cdf.a(x.star)
    b.star <- cdf.b(x.star)
    toplot <- data.frame(x = rep(x.range, 2), y = c(cdf.a(x.range), cdf.b(x.range)), d = c(rep(la, length(x.range)), rep(lb, length(x.range))))
    
    ks.plot <- ggplot(toplot) + 
    geom_step(aes(x = x, y = y, color = as.factor(d)), size = 1) +
    labs(x = 'Score Value', y = 'Cumulative Distribution') +
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    annotate('segment', x = x.star, xend = x.star, y = a.star, yend = b.star, color = 'blue', size = 1) +
    annotate('point', x = x.star, y = a.star, color = 'blue', size = 3) + 
    annotate('point', x = x.star, y = b.star, color = 'blue', size = 3)
    
    list(ks.stat = ks.stat, ks.plot = ks.plot)

}
