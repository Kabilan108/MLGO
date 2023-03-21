
# ----------------------------------- SETUP --------------------------------- #

# Set CRAN mirror
options(repos = c(CRAN = "http://lib.stat.cmu.edu/R/CRAN/"))

# Install packages
cran <- c("dplyr", "stringr", "tibble", "yaml", "doSNOW", "BiocManager")
bioc <- c("limma", "getDEE2", "GEOquery", "AnnotationDbi", "GOfuncR",
	           "org.Mm.eg.db", "biomaRt")
cran <- cran[!(cran %in% installed.packages()[, "Package"])]
bioc <- bioc[!(bioc %in% installed.packages()[, "Package"])]
if (length(cran)) install.packages(cran)
if (length(bioc)) BiocManager::install(bioc)

# Load packages
suppressMessages( suppressWarnings({
    require(dplyr)
    require(stringr)
}) )

# Load configuration
config <- yaml::read_yaml("config/config.yaml")
x <- lapply(config$path, dir.create, recursive = TRUE, showWarnings = FALSE)


# ----------------------------------- DATA ---------------------------------- #

message("Downloading metadata and annotations...")

# Download and filter the metadata
# -- For now, only use datasets with 1000 or fewer samples
metadata <- getDEE2::getDEE2Metadata(config$species$name) %>%
    filter(GEO_series != "", str_detect(Experiment_title, "^GSM")) %>%
    mutate(GSM_accession = str_extract(Experiment_title, "GSM[0-9]+")) %>%
    arrange(GEO_series) %>%
    group_by(GEO_series) %>%
    filter(n() > 1) %>%
    ungroup()

# Retrieve a biomaRt
ensembl <- biomaRt::useMart(
    "ensembl", dataset = config$species$ensembl, host = "useast.ensembl.org"
)
symbols <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl
) %>%
    rename(gene_id = ensembl_gene_id, symbol = external_gene_name) %>%
    filter(!is.na(symbol))

# Save metadata and biomaRt
save(metadata, symbols, file = paste0(config$path$raw, "/metadata.RData"))
