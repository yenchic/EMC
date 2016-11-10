# Enhanced Mode Clustering
Performing mode clustering and visualization.

- Paper reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "A comprehensive approach to mode clustering." Electronic Journal of Statistics 10.1 (2016): 210-241.
- Maintianer: yenchic@uw.edu
- `EMC.R`: the main script for running the program.
- `Example_OliveOil.R`: an example for analyzing the Olive Oil dataset.

## EMC.R
Here is a description for the main functions.

### fms
`fms = function(data, query, h, eps=1.0e-8, max.iterations=100, cut = 0.1)`
- Fast mean shift using heirachical clustering.
- Inputs:
  - data: Input data matrix.
  - query: The mesh points that you want to apply mean shift to.
  - h: Smoothing parameter.
  - max.iterations: Maximal number of iteration for mean shift.
  - eps: The tolerance. If mean shift moves less than this value, we will consider it done.
  - cut: The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
- Outputs:
  - The mode clustering result; a list consisting of
    - label: The cluster labels for query points.
    - modes: The local modes corresponding to each label.



