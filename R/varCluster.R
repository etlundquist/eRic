#' Predictor Clustering and Cluster Summary Variable Extraction
#' 
#' Perform agglomerative (hierarchical) clustering on a set of numeric predictor variables on the basis of 
#' absolute correlation distance. Create a summary variable for each cluster to use for dimension reduction.
#' 
#' @param x (data frame) data frame containing the set of numeric predictors to be clustered
#' @param y (numeric) the intended response variable (to generate x-y correlations)
#' @param corr.min (numeric) minimum correlation required to join variables into clusters
#' @param clus.summary (character) method used to extract a summary variable for each cluster (see details)
#' @param corr.method (character) which correlation coefficient to use: \code{c('pearson', 'spearman')}
#' @param corr.use (character) how to deal with missing values: \code{c('complete.obs', 'pairwise.complete.obs')}
#' @param clus.method (character) how to calculate linkages: \code{c('complete', 'single', 'centroid')}
#' 
#' @details 
#' Variables will be agglomerated into clusters on the basis of absolute correlation distance (1 - abs(cor)) 
#' using a hierarchical clustering algorithm (hclust). The resulting dendogram will be cut based on the value
#' of \code{corr.min} supplied by the user, such that all clusters contain variables with linkage correlations
#' of at least \code{corr.min}. A single summary variable will be created for each cluster based on the value
#' of \code{clus.summary} supplied by the user, with behavior detailed below. For all variables not linked into
#' a cluster, the original variable is returned as the "summary" measure unchanged.
#' \itemize{
#' \item{\strong{max.pw} - original variable with highest average pairwise correlation amongst cluster variables}
#' \item{\strong{max.y} - original variable with highest correlation with the response (\code{y})}
#' \item{\strong{avg.x} - average of standardized variables within the cluster for each observation}
#' \item{\strong{pc.1} - first principal component of cluster variables}
#' }
#' 
#' @return a list containing the following elements:
#' \itemize{
#' \item{\strong{nvar} - number of variables in \code{x}}
#' \item{\strong{nclust} - number of created clusters}
#' \item{\strong{vif} - VIF values for each original variable in \code{x}}
#' \item{\strong{clusters} - a list for each cluster containing the number of variables in the cluster, 
#' the names of the variables in the cluster, all pairwise correlations, and predictor-response correlations)}
#' \item{\strong{summaries} - a data frame containing the summary variables (one for each cluster)}
#' \item{\strong{corrplot} - a ggplot2 object visualizing the correlations amongst all predictor variables}
#' \item{\strong{params} - a list of all parameter values passed in the original call}
#' } 
#' @examples 
#' library(caret)
#' data("ChemicalManufacturingProcess")
#' x <- ChemicalManufacturingProcess[,-1]
#' y <- ChemicalManufacturingProcess[, 1]
#' res <- varCluster(x, y, corr.min = 0.75, clus.summary = 'pc.1', corr.method = 'spearman', corr.use = 'complete.obs')
#' res$nvar
#' res$nclust
#' res$clusters
#' str(res$summaries)
#' res$corrrplot
#' res$params
#' @export

varCluster <- function(x, y, corr.min = sqrt(2)/2, clus.summary = 'max.pw', corr.method = 'pearson', corr.use = 'complete.obs', clus.method = 'complete') {
    
    # load required packages
    
    if (!require(GGally) || !require(magrittr)) {
        stop('packages [GGally, magrittr] are required for this function')
    }
    
    # check input parameters
    
    if (!is.data.frame(x) || !all(sapply(x, typeof) %in% c('integer', 'double'))) {
        stop('parameter [x] must be a data frame of all numeric variables')
    }
    if (!is.numeric(y) || length(y) != nrow(x)) {
        stop('parameter [y] must be a numeric vector with length equal to nrow(x)')
    }
    if (!is.numeric(corr.min) || length(corr.min) != 1 || corr.min < 0 || corr.min > 1) {
        stop('parameter [corr.min] must be a number between [0,1]')
    }
    if (!(clus.summary %in% c('max.pw', 'max.y', 'avg.x', 'pc.1'))) {
        stop('parameter [clus.summary] must be one of [max.pw, max.y, avg.x, pc.1]')
    }
    if (!(corr.method %in% c('pearson', 'spearman'))) {
        stop('parameter [corr.method] must be one of [pearson, spearman]')
    }
    if (!(corr.use %in% c('complete.obs', 'pairwise.complete.obs'))) {
        stop('parameter [corr.use] must be one of [complete.obs, pairwise.complete.obs]')
    }
    if (!(clus.method %in% c('single', 'complete', 'centroid'))) {
        stop('parameter [clus.method] must be one of [single, complete, centroid]')
    }
    
    # calculate initial distance matrix and variable clusters
    
    cormat   <- cor(x, use = corr.use, method = corr.method)
    distmat  <- as.dist(1 - abs(cormat))
    clustree <- hclust(distmat, method = clus.method)
    clusters <- cutree(clustree, h = 1 - corr.min)
    
    # initialize root results list and store [nvar, nclust, vif] overall values

    res             <- list()
    res[['nvar']]   <- ncol(x)
    res[['nclust']] <- length(unique(clusters))
    
    efn <- function(e) {
        warning('VIF could not be calculated due to correlation matrix singularitiy')
        setNames(rep(NA, length = ncol(x)), names(x))
    }
    res[['vif']] <- tryCatch(diag(solve(cor(x, use = 'complete.obs'))), error = efn)
    
    # initialize result sub-lists
    
    res[['clusters']]  <- list()
    res[['summaries']] <- list()
    res[['params']]    <- list()
    
    # store all input parameters in a list to return
    
    res[['params']][['corr.min']]     <- corr.min
    res[['params']][['clus.summary']] <- clus.summary
    res[['params']][['corr.method']]  <- corr.method
    res[['params']][['corr.use']]     <- corr.use
    res[['params']][['clus.method']]  <- clus.method
    
    # loop through each cluster and populate cluster and summary information
    
    for (c in 1:max(clusters)) {
        
        # extract & process the set of variables mapped into each cluster
        
        cvars <- x[,names(clusters[clusters == c]), drop = FALSE]
        cname <- paste0('c', c)
        
        res[['clusters']][[cname]]              <- list()
        res[['clusters']][[cname]][['nvars']]   <- ncol(cvars)
        res[['clusters']][[cname]][['varlist']] <- names(cvars)
        
        # special handling for clusters of a single variable
        
        if (ncol(cvars) == 1) {
            res[['clusters']][[cname]][['pw.corr']] <- NA
            res[['clusters']][[cname]][['y.corr']]  <- cor(cvars[,1], y, use = corr.use, method = corr.method)
            res[['summaries']][[cname]]             <- cvars[,1]
        }
        
        # regular handling for clusters with multiple variables
        
        else {
            
            # calculate average pairwise correlations and correlations with the response
            
            pw.corr <- abs(cor(cvars, use = corr.use, method = corr.method))
            diag(pw.corr) <- NA # remove self-correlations
            pw.corr <- apply(pw.corr, 2, mean, na.rm = TRUE)
            res[['clusters']][[cname]][['pw.corr']] <- pw.corr
            
            y.corr <- apply(cvars, 2, function(x) cor(x, y, use = corr.use, method = corr.method))
            res[['clusters']][[cname]][['y.corr']] <- y.corr
            
            # calculate the cluster summary variable based on user input
            
            if (clus.summary == 'max.pw') {
                res[['summaries']][[cname]] <- cvars[,names(pw.corr[pw.corr == max(pw.corr)][1])]
            }
            else if (clus.summary == 'max.y') {
                res[['summaries']][[cname]] <- cvars[,names(y.corr[abs(y.corr) == max(abs(y.corr))][1])]
            }
            else if (clus.summary == 'avg.x') {
                cvars.std <- data.frame(lapply(cvars, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
                res[['summaries']][[cname]] <- apply(cvars.std, 1, mean, na.rm = TRUE)
            }
            else if (clus.summary == 'pc.1') {
                pcs <- prcomp(cvars[complete.cases(cvars),], retx = TRUE, center = TRUE, scale = TRUE)
                res[['summaries']][[cname]] <- rep(NA, nrow(cvars))
                res[['summaries']][[cname]][complete.cases(cvars)] <- pcs$x[,1]
            }
            else {
                stop('invalid value selected for [clus.summary]')
            }
            
        } # end if  (single-variable clusters)
    }     # end for (main cluster loop)
    
    # create corrplot object and return results
    
    toplot <- x[,names(sort(clusters))]
    names(toplot) <- abbreviate(names(sort(clusters)), minlength = 6)
    res[['corrplot']]  <- ggcorr(toplot, c(corr.use, corr.method), hjust = 1, size = 100/ncol(x), 
                                 layout.exp = 2, low = 'blue', mid = 'white', high = 'red')
    
    res[['summaries']] <- data.frame(res[['summaries']])
    res
    
}
