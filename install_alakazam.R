library(devtools)
library(versions)

d <- read.dcf("DESCRIPTION", fields = "Imports")
d <- sub("^\\n","",d)
imports <- strsplit(d,",\n")[[1]]
alakazam <- imports[grep("alakazam", imports)]

required_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", alakazam)
devel_999 <- length(grep("\\.999$",required_version)) > 0
cran_versions <- available.versions("alakazam")
cran_versions <- cran_versions$alakazam$version


if (required_version %in% cran_versions) {
    install.versions("alakazam",required_version)
} else {
    if (!devel_999) { warning(paste0(required_version," not found in CRAN. Used Bitbucket instead.")) }
    install_bitbucket("kleinstein/alakazam@default")
}