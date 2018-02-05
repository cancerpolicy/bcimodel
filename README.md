# bcimodel
Modeling package for the Breast Cancer Initiative at Fred Hutch

Contact: Jeanette Birnbaum, kurian@uw.edu

## Package structure
The package structure follows the somewhat oblique standard structure of R packages

Folder | Contains
------ | --------
inst/extdata | Original data files
data-raw | Scripts that convert data in __inst/extdata__ to formatted data in __data__
data | Formatted data that can be loaded via the package
R | Functions to support analyses
man | Help files for functions in __R__
analyses | Analyses
vignettes | Example vignette

## Install instructions
```r
devtools::install_github('cancerpolicy/bcimodel', build_vignettes=TRUE)
```

## Start here
Please see the Introduction vignette after loading the package:
```r
vignette('introduction', 'bcimodel')
```

Alternatively, [view via Git](https://github.com/cancerpolicy/bcimodel/blob/master/vignettes/introduction.Rmd) with slightly worse formatting.

