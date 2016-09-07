#' Bootstrapped Feature Importance via Filter/Model-Based Selection Methods
#' 
#' A normalized feature importance value is calculated for each potential predictor with
#' respect to all of the \code{methods} provided in the call. In addition, average importance
#' and average rank are calculated for each predictor across all methods. Right now only 
#' classification is supported, so \code{y} should be a binary numeric variable. Regression 
#' methods (pearson/spearman correlations, anova) may be supported in the future.
#'
#' @param df (data frame) data frame containing the response and all potential predictors
#' @param y (character) binary response variable name (must be a variable within \code{df})
#' @param methods (character) vector of methods to use for importance calculations (see details)
#' @param nboot (integer) number of bootstrap samples to use (or zero for no bootstrapping)
#' @param nbins (integer) number of bins to use for chi-squared / information value calculations
#' @param nplot (integer) number of variables to show on the final \code{varImp.plot}
#' @param nabin (logical) whether to include an additional bin for missing values in chi2/iv
#' @param control (list) parameters to pass to each modeling function to override/augment defaults (see details)
#' 
#' @details 
#' Using the \code{df} of predictors and \code{y} response variable supplied, feature importance
#' scores will be calculated for each of the \code{methods} supplied in the call. Supported methods 
#' and the associated variable importance metrics are described below:
#' 
#' \itemize{
#' \item{\strong{iv} - bootstrap-averaged Information Value based on binned predictor values}
#' \item{\strong{chi2} - bootstrap-averaged Chi-Squared based on binned predictor values} 
#' \item{\strong{rf} - MeanDecreaseAccuracy variable importance metric from a single RandomForest model}
#' \item{\strong{gbm} - RelativeInfluence variable importance metric from a single GBM model}
#' \item{\strong{be} - bootstrap-averaged GCV reduction variable importance metric from MARS/Earth models}
#' \item{\strong{bl} - bootstrap-averaged univariate AUC for the set of predictors selected by Lasso models}
#' }
#' 
#' As the RF/GBM methods already include inherent bootstrapping in the tree ensembles each
#' of those models is run only once. Importance scores from the other methods are derived
#' via bootstrap averaging to reduce variance and increase the stability of importance metrics.
#' If you want to specify alternate parameters for each of the modeling methods, you can pass
#' them as lists within the main \code{control} parameter. For each method you want to use, pass
#' named arguments as list items within a list where the outer name matches the method name (e.g. 'rf') 
#' See the examples for how this works. Any method parameters passed in this way will either
#' override matching defaults or be added as additional parameters to the model. Be careful 
#' with overriding the defaults as all combinations of parameters have not been fully tested!
#' 
#' @return a list containing the following elements:
#' \itemize{
#' \item{\strong{varImp.df} - a data frame containing average importance, average rank, and method-specific importance for all predictors}
#' \item{\strong{varImp.plot} - a ggplot2 object showing average normalized importance across all methods for the top \code{nplot} predictors} 
#' \item{\strong{methods} - the character vector of methods passed in the call}
#' \item{\strong{params} - a list containing the additional parameters used for each model}
#' } 
#' @examples 
#' library(caret)
#' data(GermanCredit, package = 'caret')
#' 
#' credit <- GermanCredit
#' credit$Class <- as.numeric(credit$Class == 'Good')
#' credit <- credit[,-nearZeroVar(credit)]
#' credit <- credit[,-findCorrelation(cor(select(credit, -Class)), cutoff = 0.8)]
#' 
#' res <- bootImp(credit, 'Class', nboot = 10, nbins = 10, nplot = 20)
#' res$varImp.df
#' res$varImp.plot
#' res$methods
#' res$params
#' 
#' controls <- list('rf' = list(ntree = 200, mtry = 5, nodesize = 10, importance = TRUE),
#'                 'gbm' = list(n.trees = 100, shrinkage = 0.25, cv.folds = 10)
#'                 )
#' res <- bootImp(credit, 'Class', control = controls)
#' res$varImp.df
#' res$varImp.plot
#' res$methods
#' res$params
#' @export

bootImp <- function(df, y, methods = c('iv', 'chi2', 'rf', 'gbm', 'be', 'bl'), nboot = 10, nbins = 10, nplot = 25, nabin = FALSE, control = list()) {
    
    # check input parameters
    #-----------------------
  
    if (!is.data.frame(df)) {
      stop('argument [df] must be a data frame')
    }
    if (!is.character(y) || length(y) != 1 || length(intersect(y, names(df))) == 0) {
      stop('argument [y] must the name of the response variable within the [df]')
    }
    if (length(setdiff(methods, c('iv', 'chi2', 'rf', 'gbm', 'be', 'bl'))) > 0) {
        stop('currently supported importance methods are [iv, chi2, rf, gbm, be, bl]')
    }
    if (!is.numeric(nboot) || nboot != round(nboot) || nboot < 0) {
      stop('argument [nboot] must be either a positive integer or zero (no bootstrapping)')
    }
    if (!is.numeric(nbins) || nbins != round(nbins) || nbins <= 0) {
      stop('argument [nbins] must be a positive integer (default = 10)')
    }
    if (!is.logical(nabin)) {
      stop('argument [nabin] must be logical (TRUE/FALSE)')
    }
  
    # load required packages
    #-----------------------
    
    if (!require(dplyr) || !require(ggplot2) || !require(caret)) {
        stop('packages [dplyr, ggplot2, caret] are required for this function')
    }
    if ('rf' %in% methods && !require(randomForest)) {
        stop('package [randomForest] is required for feature selection based on Random Forest (rf)')
    }
    if ('gbm' %in% methods && !require(gbm)) {
        stop('package [gbm] is required for feature selection based on GBM (gbm)')
    }
    if ('bl' %in% methods && !require(glmnet)) {
        stop('package [glmnet] is required for feature selection based on bagLasso (bl)')
    }
    
    # set up bootstrap indicators as list elements
    #---------------------------------------------
    
    nrows <- nrow(df)  
    bootind <- list()
    
    if (nboot == 0) {
        bootind[[1]] <- 1:nrows
    }
    else {
        for (b in 1:nboot) {
            bootind[[b]] <- sample(nrows, nrows, replace = TRUE)
        }
    }
    
    # initialize params & results list
    #---------------------------------
    
    params  = list()
    var.imp = list()
    
    # define a few utility subroutines shared across methods
    #-------------------------------------------------------
    
    # bin predictor using either raw values/factor levels or quantile-based bins
    bin.x <- function(x) {
        if (is.factor(x) || length(unique(x)) <= nbins) {
            bins  <- as.factor(x)
            nbins <- length(levels(bins))
        }
        else {
            bins  <- cut(x, unique(quantile(x, seq(0, 1, 1/nbins), na.rm = TRUE)), include.lowest = TRUE)
            nbins <- length(levels(bins))
        }
        if (nabin == TRUE && length(x[is.na(x)])>0) {
            levels(bins)[nbins+1] <- 'NA'
            bins[is.na(bins)]     <- 'NA'
        }
        bins
    }
    
    # produce a vector of factor-contrast expanded names (e.g. all value suffixes)
    expand.levels <- function(name, df) {
      if (is.factor(df[[name]])) paste0(name, rownames(contrasts(df[[name]])))
      else name
    }
    
    # map each factor-contrast expanded name back to the base factor name
    get.basevar <- function(var, expanded) {
      for (i in 1:length(expanded)) {
        if (length(intersect(var, expanded[[i]])) > 0) {
          return(names(expanded)[i])
        }
      }
    }
    
    # calculate Information Value (iv) if requested
    #----------------------------------------------
    
    if ('iv' %in% methods) {
    
        # calculate IV for a single x predictor
        iv.x <- function(x, y) {
            bins  <- bin.x(x)
            px_y0 <- table(bins[y==0])/length(bins[y==0])
            px_y1 <- table(bins[y==1])/length(bins[y==1])
            woe   <- ifelse(px_y0>0 & px_y1>0, log(px_y1/px_y0), NA)
            sum((px_y1 - px_y0) * woe, na.rm = TRUE)
        }
        
        # calculate IV for all predictors with a given bootstrap sample
        iv.boot <- function(b) {
            x.b <- df[b, setdiff(names(df), y)]
            y.b <- df[b, y]
            setNames(sapply(x.b, iv.x, y = y.b), colnames(x.b))
        }
        
        # calculate IV over all bootstrap samples
        print('Calculating Feature Importance via Information Value (iv)...')
        
        iv.boot.imp <- lapply(bootind, iv.boot)
        names <- attr(iv.boot.imp[[1]], 'names')
        iv.boot.imp <- cbind(names, as.data.frame(matrix(unlist(iv.boot.imp), nrow = length(names)))) %>% mutate(names = as.character(names))
        var.imp[['iv']] <- data.frame(var = iv.boot.imp$names, iv.imp = apply(iv.boot.imp[, -1, drop=FALSE], 1, mean), stringsAsFactors = FALSE)
        
    }
    
    # calculate Chi-Squared (chi2) if requested
    #------------------------------------------
    
    if ('chi2' %in% methods) {
            
        # calculate Chi-2 for a single x predictor
        chi2.x <- function(x, y) {
            bins <- bin.x(x)
            suppressWarnings(chisq.test(table(bins[!is.na(bins)], y[!is.na(bins)])))$statistic
        }
        
        # calculate Chi-2 for all predictors with a given bootstrap sample
        chi2.boot <- function(b) {
            x.b <- df[b, setdiff(names(df), y)]
            y.b <- df[b, y]
            setNames(sapply(x.b, chi2.x, y = y.b), colnames(x.b))
        }
        
        # calculate Chi-2 over all bootstrap samples
        print('Calculating Feature Importance via Chi-Squared (chi2)...')
        
        chi2.boot.imp <- lapply(bootind, chi2.boot)
        names <- attr(chi2.boot.imp[[1]], 'names')
        chi2.boot.imp <- cbind(names, as.data.frame(matrix(unlist(chi2.boot.imp), nrow = length(names)))) %>% mutate(names = as.character(names))
        var.imp[['chi2']] <- data.frame(var = chi2.boot.imp$names, chi2.imp = apply(chi2.boot.imp[, -1, drop = FALSE], 1, mean), stringsAsFactors = FALSE)
        
    }

    # calculate RandomForest (rf) if requested
    #-----------------------------------------
    
    if ('rf' %in% methods) {
      
        df.rf <- df
        df.rf[[y]] <- as.factor(df.rf[[y]]) # response needs to be factor for RF classification
        
        rf.main    <- list(formula = as.formula(paste(y, '~', '.')), data = df.rf[complete.cases(df.rf),]) 
        rf.default <- list(ntree = 500, mtry = floor(sqrt(ncol(df.rf)-1)), importance = TRUE, do.trace = 25)
        rf.params  <- c(rf.main, rf.default)
        
        if (!is.null(control[['rf']])) {
            for (i in 1:length(control[['rf']])) {
                rf.params[[names(control[['rf']])[i]]] <- control[['rf']][[i]]
            }
        }
        
        print('Calculating Feature Importance via Random Forest (rf)...')
        rf.fit <- do.call(randomForest, rf.params)
        rf.imp <- rf.fit$importance
        var.imp[['rf']] <- data.frame(var = row.names(rf.imp), rf.imp = rf.imp[,'MeanDecreaseAccuracy'], stringsAsFactors = FALSE)
        params[['rf']]  <- rf.params[-c(1,2)]
        
        
    }
    
    # calculate GBM (gbm) if requested
    #---------------------------------
    
    if ('gbm' %in% methods) {
      
        gbm.main    <- list(formula = as.formula(paste(y, '~', '.')), data = df)
        gbm.default <- list(distribution = 'bernoulli', n.trees = 250, interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.1, cv.folds = 5, class.stratify.cv = TRUE)
        gbm.params  <- c(gbm.main, gbm.default)
                         
        if (!is.null(control[['gbm']])) {
            for (i in 1:length(control[['gbm']])) {
                gbm.params[[names(control[['gbm']])[i]]] <- control[['gbm']][[i]]
            }
        }
      
        print('Calculating Feature Importance via Gradient Boosted Machine (gbm)...')
        gbm.fit <- do.call(gbm, gbm.params)
        var.imp[['gbm']] <- invisible(summary(gbm.fit, n.trees = gbm.perf(gbm.fit, plot.it = FALSE), plotit = FALSE)) %>% 
        mutate(var = as.character(var), gbm.imp = rel.inf, rel.inf = NULL)
        params[['gbm']] <- gbm.params[-c(1,2)]
        
    }

    # calculate bagEarth (be) if requested
    #-------------------------------------
    
    if ('be' %in% methods) {
      
        be.main <- list(formula = as.formula(paste(y, '~', '.')), data = df[complete.cases(df),])
        if (nboot == 0) be.default <- list(pmethod = 'backward', degree = 1)
        else            be.default <- list(pmethod = 'backward', degree = 1, B = nboot)
        be.params <- c(be.main, be.default)

        if (!is.null(control[['be']])) {
            for (i in 1:length(control[['be']])) {
                be.params[[names(control[['be']])[i]]] <- control[['be']][[i]]
            }
        }
      
        print('Calculating Feature Importance via Bagged Earth (be)...')
        if (nboot == 0) be.fit <- do.call(earth,    be.params)
        else            be.fit <- do.call(bagEarth, be.params)
        
        be.imp <- varImp(be.fit, useModel = TRUE, scale = FALSE) # varImp uses GCV reduction as metric
        be.imp <- data.frame(expanded = rownames(be.imp), be = be.imp$Overall, stringsAsFactors = FALSE)
        
        # get both factor-expanded and base names for all predictors
        impvec    <- setNames(as.list(be.imp$expanded), be.imp$expanded)
        expanded  <- setNames(lapply(as.list(names(df)), expand.levels, df), names(df))
        basevars  <- lapply(as.list(impvec), get.basevar, expanded)
        
        # fill in unexpanded factor names for those variables not showing up in importance
        for (i in 1:length(basevars)) {
            if (is.null(basevars[[i]])) {
                basevars[[i]] <- names(basevars[i])
            }
        }

        # sum importance across levels within each factor (total GCV reduction for the factor variable)
        crosswalk <- data.frame(var = unlist(basevars), expanded = names(basevars), stringsAsFactors = FALSE)
        var.imp[['be']] <- inner_join(be.imp, crosswalk, by = 'expanded') %>% dplyr::group_by(var) %>% dplyr::summarize(be.imp = sum(be)) %>% dplyr::ungroup() 
        params[['be']]  <- be.params[-c(1,2)]
        
    }
    
    # calculate bagLasso (bl) if requested
    #-------------------------------------
    
    if ('bl' %in% methods) {
        
        bl.default <- list(family = 'binomial', alpha = 1, nlambda = 100, standardize = TRUE, nfolds = 5)
        bl.params  <- bl.default
        
        if (!is.null(control[['bl']])) {
            for (i in 1:length(control[['bl']])) {
                bl.params[[names(control[['bl']])[i]]] <- control[['bl']][[i]]
            }
        }
        
        ls.boot <- function(b, bl.params) {
          
            # construct the model matrix and response vector needed for glmnet()
            x.b <- model.matrix(as.formula('~ .'), data = model.frame(df[b, setdiff(names(df), y)], na.action = na.pass))
            y.b <- df[b, y]
            
            # add the current bootstrap data set to the passed parameter list
            bl.main   <- list(x = x.b[complete.cases(x.b),], y = y.b[complete.cases(x.b)])
            bl.params <- c(bl.main, bl.params)
            
            # extract the coefficients from the cross-validated best value for lambda model
            ls.fit.b <- do.call(cv.glmnet, bl.params)
            coeff    <- as.matrix(coef(ls.fit.b$glmnet.fit, s = ls.fit.b$lambda.min))
            nzcoeff  <- coeff[coeff[,1] != 0,]
            
            # figure out which variables were/were not used in the best fit
            expanded <- setNames(lapply(as.list(names(df)), expand.levels, df), names(df))
            wasused  <- sapply(expanded, function(x) length(intersect(x, names(nzcoeff))) > 0)
            usedvars <- names(wasused[wasused == TRUE])
            notused  <- setdiff(names(wasused[wasused == FALSE]), 'surv_7y')
            
            # calculate importance as univariate AUC for all variables selected by the lasso 
            ls.imp.b <- filterVarImp(df[b, usedvars], y.b)
            ls.imp.b <- c(setNames(ls.imp.b$Overall, rownames(ls.imp.b)), setNames(numeric(length = length(notused)), notused))
            ls.imp.b[order(names(ls.imp.b))]
            
        }
            
        # average variable importance across all bootstrap samples
        print('Calculating Feature Importance via Bagged Lasso (bl)...')
        
        ls.boot.imp <- lapply(bootind, ls.boot, bl.params)  
        names <- attr(ls.boot.imp[[1]], 'names')
        ls.boot.imp <- cbind(names, as.data.frame(matrix(unlist(ls.boot.imp), nrow = length(names)))) %>% mutate(names = as.character(names))
        var.imp[['bl']] <- data.frame(var = ls.boot.imp$names, bl.imp = apply(ls.boot.imp[, -1, drop = FALSE], 1, mean), stringsAsFactors = FALSE)
        params[['bl']]  <- bl.params

    }
    
    # merge importance from all measures together
    #--------------------------------------------
    
    var.imp.df <- var.imp[[1]]
    if (length(var.imp) > 1) {
      for (r in 2:length(var.imp)) {
          var.imp.df <- inner_join(var.imp.df, var.imp[[r]], by = 'var')
      }
    }
    
    # normalize all importance values to [0-100] and calculate average importance and rank
    #-------------------------------------------------------------------------------------
    
    # calculate average importance and average rank across all used methods
    var.imp.df[,-1]     <- lapply(var.imp.df[,-1], function(x) round(100 * (x-min(x))/(max(x)-min(x)), 1))
    var.imp.df$avg.imp  <- round(apply(dplyr::select(var.imp.df, -var), 1, mean, na.rm = TRUE), 1)
    var.imp.df$avg.rank <- round(apply(sapply(dplyr::select(var.imp.df, -var, -avg.imp), function(x) row_number(desc(x))), 1, mean, na.rm = TRUE), 1)
    var.imp.df          <- dplyr::select(var.imp.df, var, avg.rank, avg.imp, everything()) %>% arrange(avg.rank)
    
    # plot the [nplot] most important variables
    var.imp.plot <- ggplot(data = var.imp.df[1:min(nplot, nrow(var.imp.df)),]) + 
    geom_bar(aes(x = reorder(as.factor(var), avg.imp), y = avg.imp), stat = 'identity') +
    coord_flip() + theme_bw() + scale_y_continuous(breaks = seq(0, 100, 10)) +
    labs(x = 'Predictor Variable', y = 'Average Normalized Importance')
    
    # return a list with full information in the DF and a summary ggplot object
    #--------------------------------------------------------------------------
    
    list(varImp.df = as.data.frame(var.imp.df), varImp.plot = var.imp.plot, methods = methods, params = params)
    
}