#!/usr/bin/env Rscript

# Install alakazam
# Checks for development version dependency and installs from Bitbucket if required.

# Import
library(devtools)
library(versions)

shazam_version <- read.dcf("DESCRIPTION", fields = "Version")
shazam_devel <- length(grep("\\.999$", shazam_version)) > 0

d <- read.dcf("DESCRIPTION", fields = "Imports")
d <- sub("^\\n", "", d)
imports <- strsplit(d, ",\n")[[1]]
alakazam <- imports[grep("alakazam", imports)]

required_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", alakazam)
alakazam_devel <- length(grep("\\.999$", required_version)) > 0
cran_versions <- available.versions("alakazam")
cran_versions <- cran_versions$alakazam$version

cat("Installing Alakazam...\n")
if (!shazam_devel & !alakazam_devel & (required_version %in% cran_versions)) {
    install.versions("alakazam", required_version)
} else {
    if (!alakazam_devel) { 
        warning(paste0(required_version," not found in CRAN. Used Bitbucket instead.")) 
    }
    install_bitbucket("kleinstein/alakazam@default")
}