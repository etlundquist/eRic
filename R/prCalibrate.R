#' Platt-Calibrated Predicted Probabilities & Calibration Plot
#' 
#' Generate Platt-Calibrated predicted probabilities for either a \strong{calibration} or 
#' \strong{validation} set as well as a \strong{calibration plot} showing actual response proportions 
#' within each of \code{nbins} predicted probability bins given class & predicted probability vectors 
#' 
#' @param r.calib (numeric) binary calibration data response vector
#' @param p.calib (numeric) numeric calibration data predicted probability vector
#' @param r.valid (numeric-optional) binary validation data response vector
#' @param p.valid (numeric-optional) numeric validation data predicted probability vector
#' @param nbins (numeric) number of bins to create for plots
#' 
#' @details 
#' Many popular machine learning algorithms produce inaccurate predicted probabilities, especially for
#' extreme observations [i.e. pr(y|x) near (0, 1)]. Examples of this include tree ensembles (RF, GBM) and 
#' SVM models, which tend to produce predicted probabilites biased towards (0.5) because of row/col bootstrap
#' sampling bias and a reliance on difficult-to-classify observations respectively. Platt (1999) proposed an
#' adjustment in these cases where the original probabilities are used as a predictor in a single-variable
#' logistic regression (optimized for log-loss) to produce more accurate adjusted predicted probabilities. 
#' This adjustment will have no effect on model AUC (sigmoid transformations are monotonic), minimal effect on 
#' standard classification accuracy, but should yield a significant improvement in log-loss/cross-entropy,
#' which is mainly used to evaluate the quality of predicted probabilities. To get a fair estimate of the effect 
#' of calibration the data should be split into 3 parts:
#' 
#' \itemize{
#' \item{\strong{model-building} - data used to train the main ML model}
#' \item{\strong{calibration} - data used to train the calibration LR classifier on out-of-sample ML predictions}
#' \item{\strong{validation} - data used to assess the effects of calibration using dual out-of-sample predictions}
#' }
#' 
#' The calibration LR model will be built using the \code{r.calib} and \code{p.calib} vectors specified
#' in the call. To get a fair estimate of the effects of calibration, \code{p.calib} should be generated
#' via out-of-sample prediction with respect to the main ML algorithm, and (\code{r.valid}, \code{p.valid})
#' should be provided to allow dual out-of-sample calibration of the original predicted probabilities. If 
#' \code{r.valid} or \code{p.valid} are not supplied, then graph/numerical metrics will be produced with
#' respect to the calibration vectors (\code{r.calib}, \code{p.calib}), but may yield biased results since 
#' these same vectors were used to fit the calibration model. See the example below for a demonstration of 
#' how this might work in a common ML context.
#' 
#' @return a list containing the following elements:
#' \itemize{
#' \item{\strong{sample} - sample to which metrics pertain: \code{c('calibration', 'validation')}}
#' \item{\strong{raw.probs} - raw unadjusted predicted probabilities passed in the call} 
#' \item{\strong{cal.probs} - Platt-calibrated predicted probabilities}
#' \item{\strong{responses} - actual responses or class labels passed in the call}
#' \item{\strong{raw.logloss} - log-loss calculated using raw predicted probabilities}
#' \item{\strong{cal.logloss} - log-loss calculated using Platt-calibrated predicted probabilities}
#' \item{\strong{cal.model} - calibration LR model object used to generate adjusted probabilities}
#' \item{\strong{cal.plot} - calibration plot showing actual class proportions within predicted probability bins}
#' \item{\strong{prob.hist} - double-histogram showing the distribution of raw and calibrated predcited probabilities}
#' } 
#' @examples 
#' library(pROC)
#' library(caret)
#' library(Information)
#' library(randomForest)
#' set.seed(1492)
#' data("train")
#' data("valid")
#' 
#' training    <- train[,setdiff(names(train), c('UNIQUE_ID', 'TREATMENT'))]
#' training    <- training[,-nearZeroVar(training)]
#' partition   <- createDataPartition(training$PURCHASE, p = 0.8, times = 1, list = FALSE)
#' building    <- training[ partition,] # use this data to train the main model
#' calibration <- training[-partition,] # use this data to calibrate predicted probabilities
#' 
#' fit     <- randomForest(as.factor(PURCHASE) ~ ., data = building, ntrees = 1000, mtry = 5)
#' p.calib <- predict(fit, calibration, type = 'prob')[,2] # predicted probabilities for calibration
#' p.valid <- predict(fit, valid,       type = 'prob')[,2] # predicted probabilities for testing
#' res     <- prCalibrate(calibration$PURCHASE, p.calib, valid$PURCHASE, p.valid, nbins = 10)
#' 
#' res$raw.logloss
#' res$cal.logloss
#' # decrease in validation log-loss
#' 
#' sum((res$raw.probs > 0.5) == res$responses)/length(res$responses)
#' sum((res$cal.probs > 0.5) == res$responses)/length(res$responses)
#' # accuracy essentially unchanged
#'
#' roc(res$responses, res$raw.probs, auc = TRUE)
#' roc(res$responses, res$cal.probs, auc = TRUE)
#' # AUC exactly the same (monotonic transformation)
#' 
#' @references 
#' \url{http://citeseer.ist.psu.edu/viewdoc/download?doi=10.1.1.41.1639&rep=rep1&type=pdf}
#' \url{http://machinelearning.wustl.edu/mlpapers/paper_files/icml2005_Niculescu-MizilC05.pdf}
#' 
#' @export

prCalibrate <- function(r.calib, p.calib, r.valid=NULL, p.valid=NULL, nbins=10) {
    
    # require external packages
    
    if (!require(ggplot2) || !require(gridExtra)) {
        stop('packages [ggplot2, gridExtra] are required for this function')
    }
    
    # input parameter validation
    
    if (!is.numeric(r.calib) || (length(unique(r.calib)) != 2) || (min(r.calib) != 0) || (max(r.calib) != 1)) {
        stop('[r.calib] must be a binary numeric vector with no missing values')
    }
    if (!is.null(r.valid) && (!is.numeric(r.valid) || (length(unique(r.valid)) != 2) || (min(r.valid) != 0) || (max(r.valid) != 1))) {
        stop('[r.valid] must be a binary numeric vector with no missing values if included')
    }
    if (!is.numeric(p.calib) || length(p.calib[is.na(p.calib)])>0 || min(p.calib)<0 || max(p.calib)>1) {
        stop('[p.calib] must be a numeric vector with values between [0,1]')
    }
    if (!is.null(p.valid) && (!is.numeric(p.valid) || length(p.valid[is.na(p.valid)])>0 || min(p.valid)<0 || max(p.valid)>1)) {
        stop('[p.valid] must be a numeric vector with values between [0,1] if included')
    }
    if (!is.numeric(nbins) || length(nbins) != 1 || nbins != round(nbins) || nbins <= 0) {
        stop('[nbins] must be an integer value greater than zero')
    }
    
    # set predictions/responses based on user input
    
    if (is.null(r.valid) || is.null(p.valid)) {
        pred   <- p.calib
        resp   <- r.calib
        sample <- 'calibration'
    }
    else {
        pred   <- p.valid
        resp   <- r.valid
        sample <- 'validation'
    }
    
    # add/subtract epsilon to predictions at (0,1) to avoid infinite log-loss
    
    pred[pred == 0] <- 1e-8
    pred[pred == 1] <- 1 - 1e-8
    
    # fit the calibration model and calculate calibrated probabilities
    
    cmodel <- glm(y ~ x, data.frame(y = r.calib, x = p.calib), family = 'binomial')
    calibrated <- predict(cmodel, data.frame(y = resp, x = pred), type = 'response')

    # calculate visualization/return measures for original probabilities (on either calibration or validation data) 

    raw.bins <- cut(pred, nbins, include.lowest = TRUE)
    raw.xval <- tapply(pred, raw.bins, mean)
    raw.yval <- tapply(resp, raw.bins, mean) 
    raw.cali <- data.frame(method = rep('Original', nbins), x = raw.xval, y = raw.yval)
    raw.logl <- (-1/length(resp)) * sum(resp*log(pred) + (1-resp)*(log(1-pred)), na.rm = TRUE)

    # calculate needed measures using transformed probabilities
    
    cal.bins <- cut(calibrated, nbins, include.lowest = TRUE)
    cal.xval <- tapply(calibrated, cal.bins, mean)
    cal.yval <- tapply(resp, cal.bins, mean)
    cal.cali <- data.frame(method = rep('Calibrated', nbins), x = cal.xval, y = cal.yval)
    cal.logl <- (-1/length(resp)) * sum(resp*log(calibrated) + (1-resp)*(log(1-calibrated)), na.rm = TRUE)

    # create the calibration and histogram plots
    
    cPlot <- ggplot(rbind(raw.cali, cal.cali)) + 
    geom_line(aes(x = x, y = y, color = as.factor(method)), size = 0.75) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2,     size = 0.50) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))                      +
    scale_color_manual(name = '', values = c('red', 'blue'))             + 
    scale_x_continuous(breaks = seq(0, 1, 0.25))                         + 
    scale_y_continuous(breaks = seq(0, 1, 0.25))                         + 
    theme_bw() + theme(legend.position = 'none')                         +
    labs(x = '', y = 'Actual Class Proportion')        
    
    hPlot <- ggplot(data.frame(method = c(rep('Original',   length(pred)), 
                                          rep('Calibrated', length(calibrated))), 
                               prob   = c(pred, calibrated)))                                     +
    geom_histogram(aes(x = prob, fill = method), 
                   bins = nbins, color = 'black', alpha = 0.50, position = 'identity')            +
    coord_cartesian(xlim = c(0, 1)) + scale_x_continuous(breaks = seq(0, 1, 0.25))                + 
    scale_fill_manual(name = '', breaks = c('Original', 'Calibrated'), values = c('blue', 'red')) +
    theme_bw() + theme(legend.position = 'bottom')                                                +
    labs(x = 'Predicted Probability', y = 'Number of Cases') 
    
    # draw the final combined plot and return results

    suppressWarnings(grid.arrange(cPlot, hPlot))
    list(sample      = sample,
         raw.probs   = pred,
         cal.probs   = calibrated,
         responses   = resp,
         raw.logloss = raw.logl, 
         cal.logloss = cal.logl, 
         cal.model   = cmodel,
         cal.Plot    = cPlot,
         prob.hist   = hPlot)
    
}
