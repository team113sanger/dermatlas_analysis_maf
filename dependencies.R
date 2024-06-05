# CONSTANTS FOR PACKAGES - user RENV syntax (see more at https://rstudio.github.io/renv/reference/install.html#ref-examples)
REQUIRED_PACKAGES <- c(
    "dplyr@1.1.4",
    "stringr@1.5.1",
    "GetoptLong@1.0.5",
    "writexl@1.5.0",
    "ggplot2@3.5.1",
    "reshape2@1.4.4",
    "plyr@1.8.9",
    "gridExtra@2.3",
    "extrafont@0.19",
    "tidyverse@2.0.0",
    "cowplot@1.1.3",
    "optparse@1.7.5",
    "tidyr@1.3.1",
    "RColorBrewer@1.1-3"
)
BIOCONDUCTOR_VERSION <- "3.16"

# UPDATE THE RENV SETTINGS
if (requireNamespace("renv", quietly = TRUE)) {
    renv::settings$bioconductor.version(BIOCONDUCTOR_VERSION)
    renv::settings$snapshot.type("all")
} else {
  message("renv is not installed. Please install renv and try again.")
  quit(status = 1)
}
