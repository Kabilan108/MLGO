

pacman::p_load(dplyr, stringr)


# Load the config file
config <- yaml::read_yaml("config/config.yaml")

# Download the metadata  and filter out the ones without GEO accessions
metadata <- getDEE2::getDEE2Metadata(config$data$species) %>%
    filter(GEO_series != "", str_detect(Experiment_title, "^GSM")) %>%
    mutate(GSM_accession = str_extract(Experiment_title, "GSM[0-9]+")) %>%
    select(GEO_series, GSM_accession, SRP_accession, SRR_accession) %>%
    arrange(GEO_series)

# Get unique GEO series
GSEs <- unique(metadata$GEO_series)
N = length(GSEs)


# Download the data
x <- pbapply::pblapply(GSEs, function(GSE) {
    print(paste("Processing ", GSE, " (", which(GSEs == GSE), "/", N, ")"))

    # Download the GSE metadata (Do not save the file
    gse_meta <- suppressMessages(
        GEOquery::getGEO(
            GSE, destdir = config$paths$temp, GSEMatrix = F, getGPL = F
        )
    )

    # Get the GSMs
    gsms <- GEOquery::GSMList(gse_meta)

    # Can `source_name_ch1` be used to split samples for DE analysis?
    sources <- sapply(gsms, function(x) x@header$source_name_ch1)
    unique_sources <- unique(sources)
    if ( length(unique_sources) >= 2 ) {
        # Split samples into 2 groups based on unique sources
        groups <- sapply(sources, function(x) which(unique_sources == x))
        check_characs <- FALSE
    } else {
        # Split samples into 2 groups based on unique characteristics
        check_characs <- TRUE
    }

    # Get the characteristics
    characs <- sapply(gsms, function(x) x@header$characteristics_ch1)


    # Delete the GSE metadata
    unlink(paste0(config$paths$temp, GSE, ".soft.gz"))
})



# Free up some memory
rm(config)
