# devtools::install_github("javh/markr")
library(markr)
library(shazam)

# Directories
pkg_path <- "."
doc_path <- "./docs"

# Build
build_mkdocs(pkg_path, doc_path = doc_path, yaml = FALSE)