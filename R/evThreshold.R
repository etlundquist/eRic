#' Classification Probability Threshold Optimization Based on Confusion Matrix Costs/Benefits
#' 
#' Calculate an optimal probability threshold for a classification problem based on 
#' the actual class distribution, the model predicted probabilities, and the costs/benefits 
#' associated with each cell in the confusion matrix.
#'
#' @param response (numeric) vector of actual class values (as a binary numeric variable)
#' @param pprob (numeric) vector of model predicted probabilities
#' @param crMatrix (numeric) matrix containing the net revenue associated with each cell in the confusion matrix (see details)
#' @param plot.points (integer) number of threshold probabilities to evaluate when making the performance metric plot (see details)
#' 
#' @details 
#' Given actual class values and model predicted probabilities one can optimize the classification threshold with respect to 
#' the costs and benefits associated with correctly identifying positive and negative cases, and making Type I and Type II errors. 
#' This is useful when overall classification accuracy doesn't align with the goals of the modeling process and you can reasonably
#' estimate the costs and benefits associated with each cell in the confusion matrix. Given this information the optimization 
#' problem then becomes one of maximizing the expected value of a new case based on: a. the predictive power of the model; b. 
#' the actual class distribution (i.e. proportion of positive cases); c. the benefits and costs associated with the 
#' new case ending up in each cell of the confusion matrix. The expected value of a new case can be expressed as:
#' 
#' \eqn{EV = pr(P) * [TPR*R(TP) - FNR*C(FN)] + pr(N) * [TNR*R(TN) - FPR*C(FP)]}
#' 
#' Where:
#' \itemize{
#' \item{\strong{pr(P)} - proportion of positive cases}
#' \item{\strong{pr(N)} - proportion of negative cases}
#' \item{\strong{TPR}   - True Positive Rate}
#' \item{\strong{R(TP)} - Revenue/Utility associated with a True Positive}
#' \item{\strong{FNR}   - False Negative Rate}
#' \item{\strong{C(FN)} - Cost/Utility associated with a False Negative}
#' \item{\strong{TNR}   - True Negative Rate}
#' \item{\strong{R(TN)} - Revenue/Utility associated with a True Negative}
#' \item{\strong{FPR}   - False Positive Rate}
#' \item{\strong{C(FP)} - Cost/Utility associated with a False Positive}
#' }
#' You need to specify the benefits/costs associated with each confusion matrix cell in \code{crMatrix}
#' where the \strong{rows} correspond to \strong{actual class values} and the \strong{columns} correspond to 
#' \strong{predicted class values}. This implies \code{crMatrix[2,2]} is the benefit associated with correctly 
#' identifying positive cases, and \code{crMatrix[1,2]} is the cost associated with mistakenly classifying a 
#' negative case as a positive one. Diagonal entries (benefits) will typically be greater than or equal to zero, 
#' and off-diagonal entries (costs) will typically be less than or equal to zero (costs expressed as negative numbers).
#' 
#' The function will return the optimal classification threshold, the unit expected value given that threshold,
#' and a ggplot2 object containing series for Sensitivity, Specificity, Accuracy, and Normalized Expected Value
#' with respect to different probability thresholds. The expected value series is normalized to [0,1] so that it
#' can be displayed on the same plot as the other metrics.
#'
#' @return a list containing the following elements:
#' \itemize{
#' \item{\strong{best.threshold} - optimal probability threshold}
#' \item{\strong{best.ev} - unit expected value given the optimal probability threshold} 
#' \item{\strong{plot.metrics} - a plot showing various performance metrics with respect to cutoff threshold}
#' }
#' @examples 
#' library(gbm)
#' library(caret)
#' data(GermanCredit, package = 'caret')
#' 
#' credit <- GermanCredit
#' credit$Class <- as.numeric(credit$Class == 'Good')
#' credit  <- credit[,-nearZeroVar(credit)]
#' gbm.fit <- gbm(Class ~ ., data = credit, n.trees = 100, shrinkage = 0.1, cv.folds = 5, distribution = 'bernoulli')
#' pprob   <- predict(gbm.fit, n.trees = gbm.perf(gbm.fit), type = 'response')
#' 
#' crMatrix <- matrix(c(0, -2, -4, 2), nrow = 2)
#' # matrix cells: [r(TN), c(FN), c(FP), r(TP)]
#' # true negatives yield no benefit, false negatives are lost potential good customers,
#' # false positives are approved bad customers, and true positives are approved good customers
#' 
#' res <- evThreshold(credit$Class, pprob, crMatrix)
#' res$plot.metrics
#' res$best.threshold
#' res$best.ev
#' @export

evThreshold <- function(response, pprob, crMatrix, plot.points = 100) {
    
    # load required functions
    #------------------------
    
    if (!require(ggplot2)) {
        stop('package [ggplot2] is required for this function')
    }
    
    # check input parameters
    #-----------------------
    
    if (!is.numeric(response) || length(unique(response)) != 2 || min(response) != 0 || max(response) != 1) {
        stop('[response] must be a binary numeric variable with no missing values')
    }
    if (!is.numeric(pprob) || length(pprob[!is.na(pprob)]) < length(pprob) || min(pprob) < 0 || max(pprob) > 1) {
        stop('[pprob] must be a numeric variable between [0,1] with no missing values')
    }
    if (!is.matrix(crMatrix) || dim(crMatrix) != c(2,2)) {
        stop('[crMatrix] must be a 2x2 matrix with net revenues for each confusion matrix cell')
    }
    if (!is.numeric(plot.points) || plot.points != round(plot.points) || plot.points <= 0) {
        stop('[plot.points] must be an integer value greater than zero')
    }
    
    # define a function to calculate UNIT EXPECTED VALUE
    #---------------------------------------------------
    
    pEV <- function(thr, response, pprob, crMatrix) {
        crMatrix[1,1] * (sum(pprob <  thr & response == 0)/sum(response == 0)) * mean(response == 0) + # TN Benefit
        crMatrix[1,2] * (sum(pprob >= thr & response == 0)/sum(response == 0)) * mean(response == 0) + # FP Cost
        crMatrix[2,1] * (sum(pprob <  thr & response == 1)/sum(response == 1)) * mean(response == 1) + # FN Cost
        crMatrix[2,2] * (sum(pprob >= thr & response == 1)/sum(response == 1)) * mean(response == 1)   # TP Benefit
    }
    
    # generate vectors of [Sensitivity, Specificity, Accuracy, Expected Value] by Threshold
    #--------------------------------------------------------------------------------------
    
    thr  <- seq(0, 1, length.out = plot.points)
    sens <- numeric(length = plot.points)
    spec <- numeric(length = plot.points)
    acc  <- numeric(length = plot.points)
    ev   <- numeric(length = plot.points)
    
    for (p in 1:plot.points) {
        sens[p] <- sum(pprob >= thr[p] & response == 1)/sum(response == 1)
        spec[p] <- sum(pprob <  thr[p] & response == 0)/sum(response == 0)
        acc[p]  <- sum(as.numeric(pprob >= thr[p]) == response)/length(response)
        ev[p]   <- pEV(thr[p], response, pprob, crMatrix)
    }
    
    # normalize the expected value vector and find the optimal expected value
    #------------------------------------------------------------------------
    
    ev     <- (ev - min(ev))/(max(ev) - min(ev))
    ev.opt <- optimize(pEV, c(0, 1), maximum = TRUE, response, pprob, crMatrix)
    
    # plot all relevant quantities and return results
    #------------------------------------------------
    
    toplot <- data.frame(threshold = rep(thr, 4),
                         label     = c(rep('Sensitivity', plot.points), 
                                       rep('Specificity', plot.points),
                                       rep('Accuracy', plot.points),
                                       rep('Normalized Expected Value', plot.points)),
                         value     = c(sens, spec, acc, ev))
    
    metrics <- ggplot(data = toplot) + geom_line(aes(x = threshold, y = value, color = label), size = 1) +
    theme_bw() + theme(legend.position = 'bottom', legend.title = element_blank()) + 
    scale_x_continuous(breaks = seq(0, 1, 0.10)) + scale_y_continuous(breaks = seq(0, 1, 0.10)) +
    labs(x = 'Probability Cutoff Threshold', y = 'Evaluation Metric')

    list(best.threshold = ev.opt$maximum, best.ev = ev.opt$objective, plot.metrics = metrics)
    
}
