[![View crameri perceptually uniform scientific colormaps on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)

![](crameri7.0.png)

# About
A simple Matlab function for [Fabio Crameri's perceptually uniform scientific colormaps](https://www.fabiocrameri.ch/colourmaps/). 

# Usage 
`crameri` without any inputs displays the options for colormaps. 

`crameri ColormapName` sets the colormap of the current axes. 

`cmap = crameri('ColormapName')` returns a 256x3 colormap. For a visual depiction of valid colormap names, type crameri. 

`cmap = crameri('-ColormapName')` a minus sign preceeding any ColormapName flips the order of the colormap. 

`cmap = crameri(...,NLevels)` specifies a number of levels in the colormap. Default value is 256. 

`cmap = crameri(...,'pivot',PivotValue)` centers a diverging colormap such that white  corresponds to a given value and maximum extents are set using current caxis limits.  If no `PivotValue` is set, 0 is assumed. 

# Citation 
Crameri, Fabio. (2021). Scientific colour maps (7.0.1). Zenodo. [https://doi.org/10.5281/zenodo.5501399](https://doi.org/10.5281/zenodo.5501399)