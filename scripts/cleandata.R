
# ----------------------------------- SETUP --------------------------------- #

# Set CRAN mirror
options(repos = c(CRAN = "http://lib.stat.cmu.edu/R/CRAN/"))

# Install packages
cran <- c("dplyr", "stringr", "tibble", "tidyr", "yaml", "feather")
cran <- cran[!(cran %in% installed.packages()[, "Package"])]
if (length(cran)) install.packages(cran)

# Load packages
require(dplyr)
require(stringr)
require(tidyr)
require(tibble)
require(feather)

# Configuration
config <- yaml::read_yaml("config/config.yaml")
x <- lapply(config$path, dir.create, recursive = TRUE, showWarnings = FALSE)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("No DEG results file (batch-*-DEG.rds) specified.")
} else {
    degfile <- args[1]
}


# ------------------------------ DATA CLEANING ------------------------------ #

# Load results
res <- readRDS(degfile)

# Clean up matrix
rownames(res) <- NULL
missing <- res[!complete.cases(res),"GSE"]
res <- as.data.frame(res[complete.cases(res),])

# Add suffix to duplicate column names
cols <- names(res)[-c(1, length(names(res)))]
dups <- cols[duplicated(cols)]
if ( length(dups) > 0 ) {
    cols[duplicated(cols)] <- paste0(cols[duplicated(cols)], "_", seq_along(dups))
}
names(res) <- c("GSE", cols, "GO_terms")

# Convert columns to numeric
res[cols] <- lapply(res[cols], as.numeric)

# Save to feather
write_feather(res, paste0(
    config$path$processed, str_extract(degfile, "batch-\\d+"), "-clean-data.feather"
))

message("Data cleaned.")
