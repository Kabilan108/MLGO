# CRAN installations
pkgs <- c("dplyr", "stringr", "tibble", "yaml", "pbapply", "BiocManager")
pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(pkgs)) install.packages(pkgs)

# Bioconductor installations
pkgs <- c("limma", "getDEE2", "GEOquery", "AnnotationDbi", "GOfuncR", "org.Mm.eg.db", "biomaRt")
pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(pkgs)) BiocManager::install(pkgs)

# Load packages
require(dplyr)
require(stringr)

# Make aliases
assays <- SummarizedExperiment::assays

# Parameters from config
species <- "mmusculus"
species_short <- "mm"

# Download the metadata
# filter out the ones without GEO accessions, those with only one sample
metadata <- getDEE2::getDEE2Metadata(species) %>%
    filter(GEO_series != "", str_detect(Experiment_title, "^GSM")) %>%
    mutate(GSM_accession = str_extract(Experiment_title, "GSM[0-9]+")) %>%
    arrange(GEO_series) %>%
    group_by(GEO_series) %>%
    filter(n() > 1) %>%
    ungroup()

# Keep only datasets with 1000 or fewer samples
metadata <- metadata %>%
    group_by(GEO_series) %>%
    filter(n() <= 1000) %>%
    ungroup()

# Get unique GEO series
GEO_series <- unique(metadata$GEO_series)
N <- length(GEO_series)

# Get the biomaRt
ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")


# Process the datasets
x <- pbapply::pblapply(GEO_series, function(GSE) {
	print(paste("Processing ", GSE, " (", which(GEO_series == GSE), "/", N, ")"))

    # Download the GSE metadata
    gse_meta <- suppressMessages( GEOquery::getGEO(
        GSE, destdir = config$paths$temp, GSEMatrix = FALSE, getGPL = FALSE
    ) )

    # Get the GSMs
    gsms <- GEOquery::GSMList(gse_meta)

    # Retrieve sample sources and characteristics
    characs <- as.data.frame( rbind(
        sapply(gsms, function(x) x@header$source_name_ch1),
        sapply(gsms, function(x) x@header$characteristics_ch1)
    ) )
    rownames(characs) <- NULL

    # Filter out terms with 'read' or 'library' or with only one unique value
    I <- !grepl("librar[^:]*:", characs[,1]) & !grepl("read[^:]*:", characs[,1])
    J <- apply(characs, 1, function(x) length(unique(x)) > 1)
    characs <- characs[I & J, ]

    if (nrow(characs) == 0) {
        # If there are no unique characteristics, return NULL
        return(NULL)
    } else if (nrow(characs) == 1) {
        # If there is only one unique characteristic, use it as the grouping variable
        characs <- unlist(characs)
        unique_characs <- unique(characs)
        groups <- sapply(characs, function(x) which(unique_characs == x))
        groups[groups != 1] <- "B"
        groups[groups == 1] <- "A"
    } else {
        # If there are multiple unique characteristics, use the one with the fewest
        # unique values
        I <- which.min( apply(characs, 1, function(x) length(unique(x))) )
        characs <- characs[I, ] %>% unlist()

        # Assign each sample to a group
        unique_characs <- unique(characs)
        groups <- sapply(characs, function(x) which(unique_characs == x))
        groups[groups != 1] <- "B"
        groups[groups == 1] <- "A"
    }

    # Get SRR accessions
    SRRs <- metadata$SRR_accession[metadata$GEO_series == GSE]
    groups <- groups[names(groups) %in% metadata$GSM_accession[metadata$GEO_series == GSE]]

    # If there are more than 600 samples, load in chunks of 600
    # (otherwise, the memory usage will be too high)
    if (length(SRRs) > 600) {
        SRR_list <- split(SRRs, ceiling(seq_along(SRRs) / 600))
    } else {
        SRR_list <- list(SRRs)
    }

    # Download the STAR gene counts
    counts <- lapply(SRR_list, function(SRR) {
        SE <- suppressMessages( suppressWarnings( getDEE2::getDEE2(
            species, SRR, metadata = metadata, counts = "GeneCounts"
        ) ) )
        assays(SE)$counts
    }) %>%
        do.call(cbind, .)

    # Perform DE analysis
    ## Create the design matrix
    groups <- factor(groups)
    design <- model.matrix(~0 + groups)
    colnames(design) <- levels(groups)
    ## Create the contrast matrix
    contrast <- limma::makeContrasts(A - B, levels = design)
    ## Fit the model
    fit <- limma::lmFit(counts, design)
    ## Perform the contrast
    cFit <- limma::contrasts.fit(fit, contrast)
    ## Perform the empirical Bayes
    eFit <- limma::eBayes(cFit)
    ## Get the results
    deg <- limma::topTable(eFit, number = Inf) %>%
        tibble::rownames_to_column("gene_id") %>%
        mutate(is_candidate = as.numeric(adj.P.Val < 0.05))

    # Update gene symbols
    symbols <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id", values = deg$gene,
        mart = ensembl
    ) %>%
        dplyr::rename(gene_id = ensembl_gene_id, symbol = external_gene_name) %>%
        dplyr::filter(!is.na(symbol))
    deg <- left_join(deg, symbols, by = "gene_id") %>%
        drop_na(symbol)

    # Gene Ontology enrichment
    res <- suppressMessages( suppressWarnings( GOfuncR::go_enrich(
        deg[c("symbol", "is_candidate")], organismDb = "org.Mm.eg.db",
        n_randsets = 100, silent = TRUE
    ) ) )$results %>%
        filter(FWER_underrep < 0.05 | FWER_overrep < 0.05)

    # Delete temporary files
    file.remove(paste0(tempdir(),"/*"), recursive = TRUE)
    file.remove(paste0("./data/.temp/*", recursive = TRUE))

    # Return the results
    list(
        x = deg$logFC,
        y = res$node_id,
        names = deg$symbol
    )
})
