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

Some users might experience issue with installing shazam if certain dependencies aren't installed correctly. One straightforward way to ensure you have everything you need is to install via Bioconductor:

```{r}
install.packages("BiocManager")
BiocManager::install("shazam")
```

Some libraries that might require manual installation (depending on your system) includes `xml2`, `rhtslib`, `igraph`, `curl`, `zlib`. 

 - If you are using a conda environment, these libraries can be installed using the following command
 ```{sh}
 conda install -c conda-forge curl \
     zlib \
     r-xml2 \
     r-igraph \
     r-rhtslib \
     xz
 ```

 - For Debian/Ubuntu, necessary system development packages can be installed using `apt-get`. The names of these equivalent libraries are `build-essential`, `libcurl4-openssl-dev`, `libssl-dev`, `liblzma-dev`, `libxml2-dev`, `xz-utils`, `libglpk-dev`, `libhts-dev`. 
 
 After this installation, xml2 and igraph can be directly installed from apt-get with library names `r-cran-xml2`, `r-cran-igraph`.

In case if your further encounter errors in installing shazam, feel free to contact us [here](https://dowser.readthedocs.io/en/stable/#contact).

