Download
-------------------------------------------------------------------------------

The latest stable release of SHazaM can be downloaded from 
[CRAN](http://cran.rstudio.com/web/packages/shazam) or 
[Bitbucket](https://bitbucket.org/kleinstein/shazam/downloads).

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install SHazaM is via CRAN:

```R
install.packages("shazam")
```

Downloaded source builds from Bitbucket may be installed in the usual way:

```R
install.packages("shazam_x.y.z.tar.gz", repos=NULL, type="source")
```

Building Development Versions
-------------------------------------------------------------------------------

To build from the [source code](http://bitbucket.org/kleinstein/shazam),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_bitbucket("kleinstein/shazam@default")
```

Note, using `install_bitbucket` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using devtools, 
roxygen and knitr:

```R
library(devtools)
install_deps()
document()
build()
install()
```
