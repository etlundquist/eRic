varCluster <- function(x, y, corr.min = sqrt(2)/2, clus.summary = 'max.pw', corr.method = 'pearson', corr.use = 'complete.obs', clus.method = 'complete') {
    
    # calculate initial distance matrix and variable clusters
    
    cormat   <- cor(x, use = corr.use, method = corr.method)
    distmat  <- as.dist(1 - abs(cormat))
    clustree <- hclust(distmat, method = clus.method)
    clusters <- cutree(clustree, h = 1 - corr.min)
    
    # initialize results list and record input parameters

    res                <- list()
    res[['nvar']]      <- ncol(x)
    res[['nclust']]    <- length(unique(clusters))
    res[['clusters']]  <- list()
    res[['summaries']] <- list()
    res[['params']]    <- list()
    
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
        res[['clusters']][[cname]] <- list()
        
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
            
            pw.corr <- cor(cvars, use = corr.use, method = corr.method)
            diag(pw.corr) <- NA # remove self-correlations
            pw.corr <- apply(pw.corr, 2, mean, na.rm = TRUE)
            res[['clusters']][[cname]][['pw.corr']] <- pw.corr
            
            y.corr <- apply(cvars, 2, function(x) cor(x, y, use = corr.use, method = corr.method))
            res[['clusters']][[cname]][['y.corr']] <- y.corr
            
            # calculate the cluster summary variable based on user input
            
            if (clus.summary == 'max.pw') {
                res[['summaries']][[cname]] <- cvars[,names(pw.corr[abs(pw.corr) == max(abs(pw.corr))][1])]
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
            
        } # end if
    }     # end for
    
    # return summary object
    
    res[['summaries']] <- data.frame(res[['summaries']])
    res
    
}
            
library(caret)
library(AppliedPredictiveModeling)
library(dplyr)
library(Information)

data("ChemicalManufacturingProcess")
chemicals <- ChemicalManufacturingProcess

x <- select(train, -TREATMENT, -PURCHASE)
y <- train$PURCHASE

corr.min     <- sqrt(2)/2
clus.summary <- 'max.pw'
corr.method  <- 'pearson'
corr.use     <- 'complete.obs'
clus.method  <- 'complete'

res <- varCluster(x, y, corr.min = 0.5, clus.summary = 'pc.1', corr.method = 'spearman', corr.use = 'pairwise.complete.obs')


            
            
            
    
    
    
        
        names(clusters[clusters == 30])
        cvars <- x[,names(clusters[clusters == 30]), drop = FALSE]
        
        avgpw <- cor(cvars, use = corr.use, method = corr.method)
        diag(avgpw) <- NA
        avgpw <- apply(avgpw, 1, mean, na.rm = TRUE)
        maxpw <- names(avgpw[abs(avgpw) == max(abs(avgpw))])
        
        y.corr <- apply(cvars, 2, function(x) cor(x, y, use = corr.use, method = corr.method))
        maxyc <- names(y.corr[abs(y.corr) == max(abs(y.corr))])
        
        pcs <- prcomp(cvars[complete.cases(cvars),], retx = TRUE, center = TRUE, scale = TRUE)
        pc1 <- rep(NA, nrow(cvars))
        pc1[complete.cases(cvars)] <- pcs$x[,1] * -1
        
        cvars.std <- data.frame(lapply(cvars, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
        cavg      <- apply(cvars.std, 1, mean, na.rm = TRUE)
        
        pcs$x[,1]    
    

sum(distmat[is.na(distmat)])
