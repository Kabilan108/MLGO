# enrichment.R
# This script performs gene ontology enrichment analysis on the differentially
# expressed genes between tumor and normal samples in the downloaded datasets.

# Required packages
pkgs <- c("limma", "yaml", "GO.db", "hgu133plus2.db", "pbapply", "dplyr",
          "tibble")
for (pkg in pkgs) {
    if (!(pkg %in% rownames(installed.packages()))) {
        install.packages(pkg)
    }
}

# Define some aliases
hgu133plus2.db <- hgu133plus2.db::hgu133plus2.db
select <- AnnotationDbi::select
keys <- AnnotationDbi::keys

# Load configuration
config <- yaml::read_yaml("config/config.yaml")$scripts$enrichment.R

# Read the list of available datasets
datasets <- readLines(config$datasets)

# Create a mapping between the probe IDs and the gene symbols
map <- select(hgu133plus2.db, keys = keys(hgu133plus2.db),
              columns = c("PROBEID", "ENTREZID"), keytype="PROBEID")


# Loop through datasets and perform enrichment analysis
res <- pbapply::pblapply(datasets, function(dataset) {
    # Load expression matrix
    exprs <- read.table(paste0(config$datadir, "/", dataset, "_exprs.tsv"),
                        header = TRUE, row.names = 1, sep = "\t") |>
        as.matrix()

    # Map probe IDs to gene symbols
    rownames(exprs)<- map$ENTREZID[match(rownames(exprs), map$PROBEID)]

    # Remove rows with missing values
    exprs <- exprs[!is.na(rownames(exprs)), ]

    # Load phenotype data
    pData <- read.table(paste0(config$datadir, "/", dataset, "_data.tsv"),
                        header = TRUE, row.names = 1, sep = "\t") |>
        as.data.frame()

    # Define groups
    group <- factor(pData$type)

    # Create the design matrix
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)

    # Create contrast matrix
    contrast <- limma::makeContrasts(tumor - normal, levels = design)

    # Fit the model
    fit <- limma::lmFit(exprs, design)

    # Perform the contrast
    cFit <- limma::contrasts.fit(fit, contrast)

    # Apply empirical Bayes moderation of standard errors.
    # This method shrinks the standard errros towards a common value, which
    # improves the stability and accuracy of the test statistics.
    eFit <- limma::eBayes(cFit)

    # Create table with results
    deg <- limma::topTable(eFit, number = Inf) |>
        tibble::as_tibble() |>
        dplyr::mutate(ID = map$PROBEID[match(ID, map$ENTREZID)])

    # Perform enrichment analysis
    go <- limma::goana(eFit, coef = 1, species = "Hs", FDR = 0.05) |>
        tibble::rownames_to_column("GO Term") |>
        tibble::as_tibble()

    print(paste0("Finished ", dataset))
})
