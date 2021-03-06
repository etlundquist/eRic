% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prCalibrate.R
\name{prCalibrate}
\alias{prCalibrate}
\title{Platt-Calibrated Predicted Probabilities & Calibration Plot}
\usage{
prCalibrate(r.calib, p.calib, r.valid = NULL, p.valid = NULL, nbins = 10)
}
\arguments{
\item{r.calib}{(numeric) binary calibration data response vector}

\item{p.calib}{(numeric) numeric calibration data predicted probability vector}

\item{r.valid}{(numeric-optional) binary validation data response vector}

\item{p.valid}{(numeric-optional) numeric validation data predicted probability vector}

\item{nbins}{(numeric) number of bins to create for plots}
}
\value{
a list containing the following elements:
\itemize{
\item{\strong{sample} - sample to which metrics pertain: \code{c('calibration', 'validation')}}
\item{\strong{raw.probs} - raw unadjusted predicted probabilities passed in the call} 
\item{\strong{cal.probs} - Platt-calibrated predicted probabilities}
\item{\strong{responses} - actual responses or class labels passed in the call}
\item{\strong{raw.logloss} - log-loss calculated using raw predicted probabilities}
\item{\strong{cal.logloss} - log-loss calculated using Platt-calibrated predicted probabilities}
\item{\strong{cal.model} - calibration LR model object used to generate adjusted probabilities}
\item{\strong{cal.plot} - calibration plot showing actual class proportions within predicted probability bins}
\item{\strong{prob.hist} - double-histogram showing the distribution of raw and calibrated predcited probabilities}
}
}
\description{
Generate Platt-Calibrated predicted probabilities for either a \strong{calibration} or 
\strong{validation} set as well as a \strong{calibration plot} showing actual response proportions 
within each of \code{nbins} predicted probability bins given class & predicted probability vectors
}
\details{
Many popular machine learning algorithms produce inaccurate predicted probabilities, especially for
extreme observations [i.e. pr(y|x) near (0, 1)]. Examples of this include tree ensembles (RF, GBM) and 
SVM models, which tend to produce predicted probabilites biased towards (0.5) because of row/col bootstrap
sampling bias and a reliance on difficult-to-classify observations respectively. Platt (1999) proposed an
adjustment in these cases where the original probabilities are used as a predictor in a single-variable
logistic regression (optimized for log-loss) to produce more accurate adjusted predicted probabilities. 
This adjustment will have no effect on model AUC (sigmoid transformations are monotonic), minimal effect on 
standard classification accuracy, but should yield a significant improvement in log-loss/cross-entropy,
which is mainly used to evaluate the quality of predicted probabilities. To get a fair estimate of the effect 
of calibration the data should be split into 3 parts:

\itemize{
\item{\strong{model-building} - data used to train the main ML model}
\item{\strong{calibration} - data used to train the calibration LR classifier on out-of-sample ML predictions}
\item{\strong{validation} - data used to assess the effects of calibration using dual out-of-sample predictions}
}

The calibration LR model will be built using the \code{r.calib} and \code{p.calib} vectors specified
in the call. To get a fair estimate of the effects of calibration, \code{p.calib} should be generated
via out-of-sample prediction with respect to the main ML algorithm, and (\code{r.valid}, \code{p.valid})
should be provided to allow dual out-of-sample calibration of the original predicted probabilities. If 
\code{r.valid} or \code{p.valid} are not supplied, then graph/numerical metrics will be produced with
respect to the calibration vectors (\code{r.calib}, \code{p.calib}), but may yield biased results since 
these same vectors were used to fit the calibration model. See the example below for a demonstration of 
how this might work in a common ML context.
}
\examples{
library(pROC)
library(caret)
library(Information)
library(randomForest)
set.seed(1492)
data("train")
data("valid")

training    <- train[,setdiff(names(train), c('UNIQUE_ID', 'TREATMENT'))]
training    <- training[,-nearZeroVar(training)]
partition   <- createDataPartition(training$PURCHASE, p = 0.8, times = 1, list = FALSE)
building    <- training[ partition,] # use this data to train the main model
calibration <- training[-partition,] # use this data to calibrate predicted probabilities

fit     <- randomForest(as.factor(PURCHASE) ~ ., data = building, ntrees = 1000, mtry = 5)
p.calib <- predict(fit, calibration, type = 'prob')[,2] # predicted probabilities for calibration
p.valid <- predict(fit, valid,       type = 'prob')[,2] # predicted probabilities for testing
res     <- prCalibrate(calibration$PURCHASE, p.calib, valid$PURCHASE, p.valid, nbins = 10)

res$raw.logloss
res$cal.logloss
# decrease in validation log-loss

sum((res$raw.probs > 0.5) == res$responses)/length(res$responses)
sum((res$cal.probs > 0.5) == res$responses)/length(res$responses)
# accuracy essentially unchanged

roc(res$responses, res$raw.probs, auc = TRUE)
roc(res$responses, res$cal.probs, auc = TRUE)
# AUC exactly the same (monotonic transformation)

}
\references{
\url{http://citeseer.ist.psu.edu/viewdoc/download?doi=10.1.1.41.1639&rep=rep1&type=pdf}
\url{http://machinelearning.wustl.edu/mlpapers/paper_files/icml2005_Niculescu-MizilC05.pdf}
}

