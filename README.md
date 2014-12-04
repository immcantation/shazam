shm
-------------------------------------------------------------------------------
December 4, 2014  
Version 0.1

Somatic hypermutation analysis package.

Build
-------------------------------------------------------------------------------
1. Clone repository
2. Update Mercurial ignore globs with (.hgignore file):
    * .*
    * *.Rproj
    * man/*.Rd
3. Use RStudio
4. Make RStudio project from existing directory
5. install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
6. Build -> Configure Build Tools
    * Check use devtools option
    * Check use roxygen option
    * Select configure roxygen options and check everything.

Requirements
-------------------------------------------------------------------------------
* R 3.0
* External R packages
    - ggplot2
    - grid
    - plyr
    - SDMTools
    - seqinr
    - zoo