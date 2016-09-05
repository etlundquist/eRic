# eRic

An R package containing some functions I wrote over the summer while working on my projects.

## Included Functions:

1. **bootImp** - calculate variable importance/perform feature selection using bootstrap resampling and a number of different methods
  - Information Value
  - Chi2
  - Random Forest
  - GBM
  - Bagged Earth/MARS 
  - Bagged Lasso

2. **evThreshold** - calculate an optimal probability threshold for classification given costs/benefits for each confusion matrix cell
  - gives you the optimal cutoff and the unit expected value
  - inputs are vector of class values, vector of predicted probabilities, and cost/benefit matrix for each confusion matrix cell
  
3. **plotKS**
  - produces a ggplot object with a KS plot for two distributions as well as the KS statistic value
  - inputs are the two distributions (e.g. defaults & non-defaults) as well as class labels for the graph
  
4. **plotYXbin**
  - produces a ggplot object with bin values of Y with respect to bins of X
  - Y metrics include: [proportions, log-odds, WOE]
  - X split methods include: [quantile-based, uniform splits, rpart-based splits]
  - Can specify the desired number of bins and whether a missing value bin should be added
  - Additionally calculates Information Value (IV) and Chi2 Statistic for the XY relationship
  
5. **tsSummary**
  - calculate summary statistics (using passed summary functions) with respect to wide-format time series variables
  
6. **tsTrends**
  - calculate linear trends with respect to wide-format time series variables

