# eRic

An R package containing functions I've written while working on predictive modeling projects. The functions are briefly described below, but documented in more detail in the function-specific help (.rd) files. All content can be downloaded/installed and then used like any other package after following the installation instructions below.

## Installation Options

Install via `devtools` straight from GitHub:

```
install.packages("devtools")
library("devtools")
install_github("etlundquist/eRic")
```

Download the tarball `eRic_0.0.0.9000.tar.gz` and install from source:

```
install.packages('/path/to/tarball/eRic_0.0.0.9000.tar.gz', repos = NULL)
```

## Included Functions

1. **bootImp** - calculate variable importance/perform feature selection using bootstrap resampling and a combination of different filter/model-based methods. Available methods are:
  - Information Value
  - Chi2
  - Random Forest
  - GBM
  - Bagged Earth/MARS 
  - Bagged Lasso

2. **evThreshold** - calculate an optimal probability threshold for classification given the costs/benefits for each confusion matrix cell
  - provides an estimate of the optimal probability cutoff with respect to confusion matrix utility
  - provides an estimate of unit expected value given your model and the optimal probability cutoff
  
3. **varCluster** - cluster predictor variables and extract cluster summaries for dimension reduction
  - use agglomerative clustering (hclust) to group highly correlated sets of predictor variables
  - create a single variable to summarize each cluster (highest pairwise, highest x-y, centroid, PC1)
  - produce a correlation plot to visualize correlation structure in predictor matrix

4. **plotYXbin** - produces a ggplot object with bin values of Y with respect to bins of X
  - Y metrics include: [proportions, log-odds, WOE]
  - X split methods include: [quantile-based, uniform splits, rpart-based splits]
  - Can specify the desired number of bins and whether a missing value bin should be added
  - Additionally calculates Information Value (IV) and Chi2 Statistic for the XY relationship
  
5. **prCalibrate** - perform Platt Scaling on raw model predicted probabilities to better align with actual class proportions and produce a calibration plot to visualize results
  - scale predicted probabilities with respect to either a calibration or independent validation set
  - produce a calibration plot showing the relationship between actual and predicted class proportions
  
6. **plotKS** - produces a ggplot object with a KS plot for two distributions as well as the KS statistic value
 
7. **tsSummary** - calculate summary statistics (using passed summary functions) with respect to wide-format time series variables
  
8. **tsTrends** - calculate linear trends with respect to wide-format time series variables

