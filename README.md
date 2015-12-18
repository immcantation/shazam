shm
-------------------------------------------------------------------------------
December 18, 2015  
Version 0.1.1

Somatic hypermutation analysis package.

Dependencies
-------------------------------------------------------------------------------
R 3.0  
R packages

  - alakazam
  - data.table
  - doParallel
  - dplyr
  - foreach
  - ggplot2
  - iterators
  - scales  
  - SDMTools
  - seqinr
  - tidyr

Build Instructions
-------------------------------------------------------------------------------
Install build dependencies:
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

Building with Rstudio:

-  Build -> Configure Build Tools
-  Check use devtools option
-  Check use roxygen option
-  Select configure roxygen options and check everything.
-  Build -> Build and Reload

Building from the R console:

```R
library(roxygen2)
library(devtools)
install_deps()
document()
build(vignettes=FALSE)
install()
```