# Download & Installation

Download
-------------------------------------------------------------------------------

The latest stable release of SHazaM can be downloaded from <a href="http://cran.rstudio.com/web/packages/shazam" target="_blank">CRAN</a>
 or <a href="https://github.comm/immcantation/shazam/tags" target="_blank">GitHub</a>.

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install SHazaM is via CRAN:

```R
install.packages("shazam")
```

Downloaded source builds from GitHub may be installed in the usual way:

```R
install.packages("shazam_x.y.z.tar.gz", repos=NULL, type="source")
```

Building Development Versions
-------------------------------------------------------------------------------

To build from the [source code](http://github.com/immcantation/shazam),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_github("immcantation/shazam@master")
```

Note, using `install_github` will not build the documentation. To generate the
documentation, clone the repository and build as normal using devtools,
roxygen and knitr:

```R
library(devtools)
install_deps()
document()
build()
install()
```

Some users might experience issues with shazam if certain dependencies are not installed correctly. One straightforward way to ensure all required packages are available is to install install shazam via Bioconductor:

```{r}
install.packages("BiocManager")
BiocManager::install("shazam")
```

If you still encouner issues with missing dependencies, you can use the [Posit Public Package Manager](https://packagemanager.posit.co/client/#/) to identify and install them manually. 


If problems persist while installing shazam, feel free to contact us [here](https://immcantation.readthedocs.io/en/stable/about.html) for support.
