#!/bin/bash

# Update package list and install essential build tools
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev

# Install R and required packages
sudo apt-get install -y r-base r-base-dev
sudo R -e "install.packages(c('dplyr', 'stringr', 'tibble', 'yaml', 'doSNOW', 'feather', 'BiocManager'), repos='http://lib.stat.cmu.edu/R/CRAN/')"
sudo R -e "BiocManager::install(c('limma', 'getDEE2', 'GEOquery', 'AnnotationDbi', 'GOfuncR', 'org.Mm.eg.db', 'biomaRt'))"

# Install Python and required packages
sudo apt-get install -y python3 python3-pip python3-venv
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install scikit-learn pandas numpy feather-format tqdm PyYAML

# Install additional system dependencies
sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev

# Set up data directories (adjust paths as needed)
mkdir -p \
    /home/ubuntu/mlgo-data/raw \
    /home/ubuntu/mlgo-data/processed \
    /home/ubuntu/mlgo-data/final \
    /home/ubuntu/mlgo-data/.temp \
    /home/ubuntu/mlgo-data/logs

# Clone repository and navigate to directory (adjust URL and directory name as needed)
git clone https://github.com/Kabilan108/MLGO.git
cd MLGO
