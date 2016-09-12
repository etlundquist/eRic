prCalibrate <- function(responses, predictions, nbins=10, cmodel=NULL) {
    
    # require external packages
    
    if (!require(ggplot2) || !require(gridExtra)) {
        stop('packages [ggplot2, gridExtra] are required for this function')
    }
    
    # parameter validation
    
    if (!is.numeric(responses) || (length(unique(responses)) != 2) || (min(responses) != 0) || (max(responses) != 1)) {
        stop('[responses] must be a binary numeric vector with no missing values')
    }
    if (!is.numeric(predictions) || length(predictions[is.na(predictions)])>0 || min(predictions)<0 || max(predictions)>1) {
        stop('[predictions] must be a numeric vector with values between [0,1]')
    }
    if (!is.numeric(nbins) || length(nbins) != 1 || nbins != round(nbins) || nbins <= 0) {
        stop('[nbins] must be an integer value')
    }
    if (!is.null(cmodel) && (length(intersect('glm', class(cmodel))) == 0 || cmodel$family$family != 'binomial' || cmodel$formula != 'responses ~ predictions')) {
        stop('[cmodel] must be a binomial glm model object fit with the formula [responses ~ predictions]')
    }
    
    # calculate needed measures using original probabilities

    old.bins <- cut(predictions, nbins, include.lowest = TRUE)
    old.xval <- tapply(predictions, old.bins, mean)
    old.yval <- tapply(responses, old.bins, mean) 
    old.cali <- data.frame(method = rep('Original', nbins), x = old.xval, y = old.yval)
    old.logl <- (-1/length(responses)) * sum(responses*log(predictions) + (1-responses)*(log(1-predictions)))
    
    # calculate transformed predicted probabilities
    
    if (!is.null(cmodel)) {
        calibrated <- predict(cmodel, data.frame(responses = responses, predictions = predictions), type = 'response')
    }
    else {
        calibrated <- glm(responses ~ predictions, family = 'binomial')$fitted.values
    }
    
    # calculate needed measures using transformed probabilities
    
    cal.bins <- cut(calibrated, nbins, include.lowest = TRUE)
    cal.xval <- tapply(calibrated, cal.bins, mean)
    cal.yval <- tapply(responses, cal.bins, mean)
    cal.cali <- data.frame(method = rep('Calibrated', nbins), x = cal.xval, y = cal.yval)
    cal.logl <- (-1/length(responses)) * sum(responses*log(calibrated) + (1-responses)*(log(1-calibrated)))

    # create the calibration and histogram plots
    
    cPlot <- ggplot(rbind(old.cali, cal.cali)) + 
    geom_line(aes(x = x, y = y, color = as.factor(method)), size = 0.75) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2,     size = 0.50) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))                      +
    scale_color_manual(name = '', values = c('red', 'blue'))             + 
    scale_x_continuous(breaks = seq(0, 1, 0.25))                         + 
    scale_y_continuous(breaks = seq(0, 1, 0.25))                         + 
    theme_bw() + theme(legend.position = 'none')                         +
    labs(x = '', y = 'Actual Class Proportion') 
    
    hPlot <- ggplot(data.frame(method = c(rep('Original',   length(predictions)), 
                                          rep('Calibrated', length(calibrated))), 
                               prob   = c(predictions, calibrated)))                              +
    geom_histogram(aes(x = prob, fill = method), 
                   bins = nbins, color = 'black', alpha = 0.50, position = 'identity')            +
    coord_cartesian(xlim = c(0, 1)) + scale_x_continuous(breaks = seq(0, 1, 0.25))                + 
    scale_fill_manual(name = '', breaks = c('Original', 'Calibrated'), values = c('blue', 'red')) +
    theme_bw() + theme(legend.position = 'bottom')                                                +
    labs(x = 'Predicted Probability', y = 'Number of Cases') 
    
    # draw the final combined plot and return results

    suppressWarnings(grid.arrange(cPlot, hPlot))
    list(old.logloss = old.logl, 
         cal.logloss = cal.logl, 
         calibrated  = calibrated, 
         cal.Plot    = cPlot,
         prob.hist   = hPlot)
    
}

library(gbm)
library(caret)
library(pROC)

data("GermanCredit")
credit <- GermanCredit
credit$Class <- as.numeric(credit$Class == 'Good')
credit <- credit[,-nearZeroVar(credit)]

set.seed(1492)
tind <- createDataPartition(credit$Class, p = 0.75, times = 1, list = FALSE)
trdf <- credit[ tind,]
vldf <- credit[-tind,]

fit   <- gbm(Class ~ ., data = trdf, n.trees = 100, shrinkage = 0.1)
probs <- predict(fit, vldf, n.trees = gbm.perf(fit), type = 'response')
res   <- prCalibrate(vldf$Class, probs, nbins = 10)

sum((probs > 0.5) == credit$Class)/length(credit$Class)
sum((res$calibrated > 0.5) == credit$Class)/length(credit$Class)

roc(credit$Class, probs, auc = TRUE)
roc(credit$Class, res$calibrated, auc = TRUE)

res$old.logloss
res$cal.logloss
                      
