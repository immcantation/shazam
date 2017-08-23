#!/usr/bin/env Rscript

# Install dependencies from CRAN or Bitbucket as needed.

# Imports
library(devtools)
library(versions)

# Function to install packages
installDep <- function(this_pack_v, dep_pack_name, dep_pack_v) {
    required_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", dep_pack_v)
    devel <- length(grep("\\.999$", required_version)) > 0
    
    cran_versions <- available.versions(dep_pack_name)
    cran_versions <- cran_versions[[dep_pack_name]]$version
    
    this_pack_devel <- length(grep("\\.999$", this_pack_v)) > 0
    
    if (!this_pack_devel & !devel & (required_version %in% cran_versions)) {
        tryCatch({ devtools::install_version(dep_pack_name, required_version, repos="https://cran.cnr.berkeley.edu") },
                 error=function(e) { 
                     cat(e, "\n")
                     message("Installing from Bitbucket...\n ")
                     install_bitbucket(paste0("kleinstein/", dep_pack_name, "@default"))
                 })
    } else {
        if (!devel) { 
            warning(paste0(required_version," not found in CRAN. Install most recent version from Bitbucket instead.")) 
        }
        install_bitbucket(paste0("kleinstein/", dep_pack_name, "@default"))
    }
}

# Parse Imports field in DESCRIPTION
pkg_version <- read.dcf("DESCRIPTION", fields="Version")
d <- read.dcf("DESCRIPTION", fields="Imports")
d <- sub("^\\n", "", d)
imports <- strsplit(d, ",\n")[[1]]

# Install alakazam
alakazam_version <- imports[grep("alakazam", imports)]
installDep(pkg_version, "alakazam", alakazam_version)