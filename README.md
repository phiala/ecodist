**ecodist R package**

Dissimilarity-based analysis functions including ordination and Mantel test functions, intended for use with spatial and community data.


**CHANGES in ecodist 2.0.10**

The spdep package overwrites the base behavior of dim() for dist objects, which breaks ecodist::MRM(). A fix has been added to ecodist 2.0.10 that returns dim.dist() to its base state, as long as ecodist is loaded after spdep. Note that this may in turn break spdep functions. The maintainer of spdep plans to remove the problematic dim.dist() eventually, since ecodist is not the only package it breaaks.

Until then, if you need ecodist and spdep both, load ecodist last. Since dim() is called by functions outside of ecodist, I can't simply specify ecodist::dim() in all relevant cases.

dim(dist(matrix(1:15, ncol=3))) 

should return NULL and ***not*** c(5, 5) if the correct (base) dim() is being used.


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

