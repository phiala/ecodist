**ecodist R package**

Dissimilarity-based analysis functions including ordination and Mantel test functions, intended for use with spatial and community data.

**CHANGES in ecodist 2.0**

  - fixed bug in crosstab() that affected expansion of single-row or -column tables using allrows or allcols; changed result to data frame
  - added icov argument to distance() for use with Mahalanobis distance
  - changed stress calculation in nmds() to match vegan and MASS calculations; formerly was a similar method that was monotonically related, but not identical.
  - added plot.nmds to display stress and r2 for NMDS ordinations across a range of dimensions
  - added addord method to add new data to an existing NMDS ordination. 
  - added clusterlevel to calculate Mantel tests for specified groupings
  - added logistic regression to MRM
  - added xdistance cross-distance function, and cross-dissimilarity analysis functions xmantel, xmgram
  - updated examples in help files to be more helpful
  - added a vignette listing dissimilarity-based analyses

