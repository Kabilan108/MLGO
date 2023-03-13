
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
    require(doSNOW)
    require(parallel)
    require(progress)
}) )

# Define aliases for functions
assays <- SummarizedExperiment::assays

# Load configuration
config <- yaml::read_yaml("config/config.yaml")
x <- lapply(config$path, dir.create, recursive = TRUE, showWarnings = FALSE)

# Get list of GSEs to download
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("No GSEs provided")
} else {
    GSE_file <- args[1]
    GEO_series <- read.table(GSE_file, header = FALSE, stringsAsFactors = FALSE)[,1]
    N <- length(GEO_series)
}

# Load metadata
if ( file.exists(paste0(config$path$raw, "/metadata.RData")) ) {
    load( paste0(config$path$raw, "/metadata.RData") )
} else {
    stop("Could not find data/processed/metadata.RData. Please run loadmetadata.R first.")
}


# ------------------------------ PARALLELIZATION ---------------------------- #

logfile <- paste0(config$path$logs, "/", str_extract(GSE_file, 'batch-\\d'), "loaddatasets.log")

CL <- makeCluster(2, outfile = logfile)
registerDoSNOW(CL)

PB <- progress_bar$new(
    format = "Processing :gse :current/:total [:bar] :percent eta: :eta",
    total = N,
    clear = FALSE,
    width = 60
)
progress <- function(i) {
    PB$tick(tokens = list(current = i, total = N, gse = GEO_series[i]))
}
opts <- list(progress = progress)


# ----------------------------------- MAIN ---------------------------------- #
start_time <- Sys.time()
res <- foreach(
    i = seq_along(GEO_series),
    .combine = rbind,
    .options.snow = opts,
    .packages = c("dplyr", "stringr", "tibble", "getDEE2", "GEOquery", "limma",
                 "GOfuncR")
) %dopar% {

    GSE <- GEO_series[i]
    message("Processing ", GSE, " (", i, "/", N, ")")

    tryCatch({
        # Download the GSE metadata
        gse_meta <- suppressMessages( GEOquery::getGEO(
            GSE, destdir = config$path$temp, GSEMatrix = FALSE, getGPL = FALSE
        ) )
        gsms <- GEOquery::GSMList(gse_meta)

        # Retrieve sample sources and characteristics
        characs <- as.data.frame( rbind(
            sapply(gsms, function(x) x@header$source_name_ch1),
            sapply(gsms, function(x) x@header$characteristics_ch1)
        ) )
        characs <- characs[
            !grepl("librar[^:]*:", characs[,1]) &
            !grepl("read[^:]*:", characs[,1]),
        ]
        characs <- characs[
            apply(characs, 1, function(x) length(unique(x)) > 1),
        ]

        if (nrow(characs) == 0) {  # No characteristics
            SKIP <- TRUE
        } else if (nrow(characs) == 1) {  # Single characteristic
            characs <- unlist(characs)
            unique_characs <- unique(characs)
            groups <- sapply(characs, function(x) which(unique_characs == x))
            groups[groups != 1] <- "B"
            groups[groups == 1] <- "A"
            SKIP <- FALSE
        } else {  # Multiple characteristics
            characs <- unlist( characs[
                which.min( apply(characs, 1, function(x) length(unique(x))) ),
            ] )
            unique_characs <- unique(characs)
            groups <- sapply(characs, function(x) which(unique_characs == x))
            groups[groups != 1] <- "B"
            groups[groups == 1] <- "A"
            SKIP <- FALSE
        }

        if (SKIP) {
	    print('I AM SKIPPING THIS DATASET')
            c(GSE, NA, NA)
        } else {
            # Get SRR accessions
            SRRs <- metadata$SRR_accession[metadata$GEO_series == GSE]
            groups <- groups[
                names(groups) %in% metadata$GSM_accession[metadata$GEO_series == GSE]
            ]

            # Split SRRs into chunks of 600
            if (length(SRRs) > 600) {
                SRR_list <- split(SRRs, ceiling(seq_along(SRRs) / 600))
            } else {
                SRR_list <- list(SRRs)
            }

            # Download counts
            counts <- lapply(SRR_list, function(SRR) {
                SE <- suppressMessages( suppressWarnings( getDEE2::getDEE2(
                    config$species$name, SRR, metadata = metadata,
                    counts = "GeneCounts", quiet = TRUE
                ) ) )
                assays(SE)$counts
            }) %>%
                do.call(cbind, .)

	    # Delete downloaded files
            unlink(paste0(config$path$temp, "/", GSE, "/*"))

            # Perform DE analysis
            groups <- factor(groups)
            design <- model.matrix(~0 + groups)
            colnames(design) <- levels(groups)
            contrast <- limma::makeContrasts(A - B, levels = design)
            fit <- limma::lmFit(counts, design)
            fit <- limma::contrasts.fit(fit, contrast)
            fit <- limma::eBayes(fit)
            deg <- limma::topTable(fit, number = Inf) %>%
                tibble::rownames_to_column("gene_id") %>%
                mutate(is_candidate = as.numeric(adj.P.Val < 0.05))

            if ( !any(deg$is_candidate == 1) ) {
		print("NO CANDIDATES")
                c(GSE, NA, NA)
            } else {
                # Update gene symbols
                deg <- left_join(deg, symbols, by = "gene_id") %>%
                    tidyr::drop_na(symbol)

                # Gene Ontology enrichment
                res <- suppressMessages( suppressWarnings( GOfuncR::go_enrich(
                    deg[c("symbol", "is_candidate")], organismDb = "org.Mm.eg.db",
                    n_randsets = 100, silent = TRUE
                ) ) )$results %>%
                    filter(FWER_underrep < 0.05 | FWER_overrep < 0.05)

                # Store results in row vector
                row <- c(GSE, deg$logFC, paste(res$node_id, collapse = ";"))
                names(row) <- c("GSE", deg$symbol, "GO_terms")
                row
            }
        }
    }, error = function(e) {
	print("SOME OTHER ERROR")
        c(GSE, NA, NA)
    })
}
end_time <- Sys.time()
message("\nElapsed time: ", difftime(end_time, start_time, units="auto"))

# --------------------------------- CLEAN UP -------------------------------- #

stopCluster(CL)
saveRDS(res, paste0(
    config$path$processed, "/", str_extract(GSE_file, 'batch-\\d'), "-DEG.rds"
))
