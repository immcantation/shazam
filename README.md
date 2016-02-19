shazam
-------------------------------------------------------------------------------
February 19, 2016
Version 0.1.2

Somatic hypermutation analysis package.

Dependencies
-------------------------------------------------------------------------------
R 3.1.2  
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
  - stringi
  - tidyr

Build Instructions
-------------------------------------------------------------------------------
Install build dependencies:
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

Building with Rstudio:

- _Build_ -> _Configure Build Tools_
- Check the _Use devtools package functions_ option
- Check the _Generate documentation with Roxygen_ option
- Select _Configure..._ Roxygen options and check everything.
- _Build_ -> _Build and Reload_

Building from the R console:

```R
devtools::install_deps()
devtools::document()
devtools::build()
devtools::install()
```

Optionally, you can skip the vignettes:
```R
devtools::build(vignettes=FALSE)
```
