#!/bin/bash
#
# This script will download over 11000 datasets from the DEE2 repository using the
# `loaddatasets.R` script.
# Datasets will be loaded in batches of ~3500 datasets. After each batch, the
# `cleandata.R` script will be run to clean the data format the results as dataframes.
# Once all 3 batches are complete, the `prepdata.py` script will be run to prepare
# the feature matrix and target vector for the machine learning models.

for file in data/raw/batch-*-GSE.txt; do
    batch=${file#data/raw/batch-}
    batch=${batch%-GSE.txt}

    echo "Processing Batch $batch"

    Rscript loaddatasets.R "$file"

    if [[ -e data/processed/DEG_B${batch}.rds ]]; then
        Rscript cleandata.R "data/processed/DEG_B${batch}.rds"
    else
        echo "Error: Batch $batch failed"
        echo "data/processed/DEG_B${batch}.rds not found for batch $batch"
    fi
done
