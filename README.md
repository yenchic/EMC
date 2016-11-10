# Enhanced Mode Clustering
Performing mode clustering and visualization.

- Paper reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "A comprehensive approach to mode clustering." Electronic Journal of Statistics 10.1 (2016): 210-241.
- Maintianer: yenchic@uw.edu
- `EMC.R`: the main script for running the program.
- `Example_OliveOil.R`: an example for applying `EMC` method to the Olive Oil dataset.

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

### EMC
`EMC = function(data, h=NULL, eps=1.0e-8, max.iterations=100, ...) UseMethod("EMC")`
`EMC.default = function(data, h=NULL, eps=1.0e-8, max.iterations=100, n0= NULL, rho=NULL, cut=0.1, noisy = F, T_denoise =5)`
- Enhanced mode clustering.
- Inputs:
  - data: Input data matrix.
  - h: Smoothing parameter.
  - max.iterations: Maximal number of iteration for mean shift.
  - eps: The tolerance. If mean shift moves less than this value, we will consider it done.
  - n0: The thresholding size for tiny clusters. Default is to use the method given in Chen et al. (2014).
  - rho: The contrast parameter for visualization. Default is to use the method given in Chen et al. (2014).
  - noisy: True or False. To desplay noisy clusters (without thresholding). Default is False.
  - T_denoise: Maximal number of denoising. If tiny clusters presence, we will remove them and redo mode clustering. This is the maximal number of redoing mean shift clustering.
  - cut: The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
- Output:
  - An S4 object about summary informations using enhanced mode clustering. A list consisting:
    - label: The cluster labels for query points.
    - modes: The local modes corresponding to each label.
    - c.matrix: The connectivity matrix.
    - vis.data: The visualization coordinates for data points.
    - vis.modes: The visualization coordinates for local modes.
    - SC.plot: The size of ordered clusters before denoising.
    - size.threshold: The size threshold for denoising tiny clusters.
    - bandwidth: The smoothing bandwidth.
    - rho: The contrast paramter used for visualization.
    - noisy.label: The cluster labels for query points before denoising.
    - noisy.modes: The local modes corresponding to each label before denoising.


