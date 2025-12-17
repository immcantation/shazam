#!/usr/bin/env Rscript

# Install dependencies from CRAN or GitHub as needed.

####
# Add here Biocondutor dependencies that
# are not installed in the immcantation/test container
bioconductor_deps <- NULL
####

library(devtools)
library(versions)
# All immcantation R packages
immcantation_packages <- c("alakazam", "shazam", "tigger", "scoper", "dowser", "enchantr")
# 
document <- list(
    "alakazam" = TRUE,
    "shazam" = FALSE,
    "tigger" = FALSE,
    "scoper" = FALSE,
    "dowser" = FALSE,
    "enchantr" = FALSE
)

# Function to install dependencies
installDep <- function(pkg, devel_mode, immcantation = immcantation_packages,
                       repos = "http://lib.stat.cmu.edu/R/CRAN/") {
    # Required version
    pkg_name <- strsplit(pkg, " ")[[1]][1]
    pkg_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", pkg)
    pkg_logic <- gsub("(.*)( *)\\(([><=]*)[^0-9.]*(.*)\\)$", "\\3", pkg)
    if (pkg_version == pkg_name) {
        pkg_version <- NULL
        pkg_logic <- ">="
    }

    is_immcantation <- pkg_name %in% immcantation

    message("Package requested: ", pkg)
    message("devel_mode: ", devel_mode)
    message("is_immcantation: ", is_immcantation)

    # Installed packages
    is_installed <- F
    installed_version <- NA
    installed_packages <- installed.packages()
    installed_pkg <- installed_packages[installed_packages[, "Package"] == pkg_name, , drop = F]

    if (nrow(installed_pkg) > 0) {
        # Check version of installed package
        tmp_pkg_version <- pkg_version
        if (is.null(pkg_version)) {
            tmp_pkg_version <- "0.0.0"
        }
        installed_version <- numeric_version(installed_pkg[, "Version"])
        is_installed <- eval(parse(text = paste0(
            "numeric_version('", installed_version, "') ",
            pkg_logic,
            " numeric_version('", tmp_pkg_version, "')"
        )))
    }

    if (!is_immcantation | !devel_mode) {
        if (is_installed) {
            message(pkg, " is available.")
        } else {
            # Install from CRAN
            tryCatch(
                {
                    devtools::install_version(pkg_name, paste(pkg_logic, pkg_version), repos = "http://lib.stat.cmu.edu/R/CRAN/", upgrade = "never")
                },
                error = function(e) {
                    # This is needed if there is an Immcantation release package that is not
                    # available from CRAN
                    cat(as.character(e), "\n")
                    message("Installing ", pkg, " from GitHub...\n ")
                    install_github(paste0("immcantation/", pkg_name, "@", pkg_version), build = TRUE)
                }
            )
        }
    } else {
        if (!document[[pkg_name]]) {
            message(paste0(pkg, ": installing most recent version from GitHub @master."))
            install_github(paste0("immcantation/", pkg_name, "@master"), upgrade = "never", force = FALSE)
        } else {
            message(paste0(pkg, ": installing and documenting most recent version from GitHub @master."))
            pkg_tmp_dir <- file.path(tempdir(), pkg_name)
            pkg_repo <- paste0("immcantation/", pkg_name, "@master")
            dir.create(pkg_tmp_dir, recursive = TRUE, showWarnings = FALSE)
            # Clone the repository
            tryCatch(
                {
                    # Clone the repository
                    system2("git", c(
                        "clone", "--branch", "master", "--depth", "1",
                        paste0("https://github.com/immcantation/", pkg_name, ".git"),
                        pkg_tmp_dir
                    ))

                    # Install dependencies, document, build and install
                    devtools::install_deps(pkg_tmp_dir, dependencies = TRUE, upgrade = "never", force = FALSE)
                    devtools::document(pkg_tmp_dir)
                    devtools::build(pkg_tmp_dir)
                    devtools::install(pkg_tmp_dir, upgrade = "never", force = FALSE)
                    unlink(pkg_tmp_dir, recursive = TRUE)
                    message(paste0("Successfully installed ", pkg, " from GitHub master"))
                },
                error = function(e) {
                    message(paste0("Error installing and documenting", pkg, " from GitHub: ", e$message))
                }
            )
        }
    }
}

# Parse this package version in DESCRIPTION
this_pkg_version <- read.dcf("DESCRIPTION", fields = "Version")
this_pkg_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", this_pkg_version)

# If the package is using the devel version number (ends in .999)
# always install devel versions of immcantation packages
devel_mode <- grepl("\\.999$", this_pkg_version)

# Parse Imports field in DESCRIPTION
d <- read.dcf("DESCRIPTION", fields = "Imports")
d <- sub("^\\n", "", d)
imports <- strsplit(d, ",\n")[[1]]

# Install immcantation packages first
immcantation <- unlist(sapply(immcantation_packages, grep, imports))
imports <- imports[unique(c(immcantation, 1:length(imports)))]
# Skip bioconductor packages
imports <- imports[!imports %in% bioconductor_deps]

# Install Bioconductor dependencies first
if (!is.null(bioconductor_deps)) {
    BiocManager::install(bioconductor_deps, update=FALSE, ask=FALSE, force=FALSE)
}

# Install
for (i in 1:length(imports)) {
    this_import <- imports[i]
    installDep(this_import, devel_mode)
}