# **MLGO:** Machine Learning for Gene Ontology Terms

## Gene Ontology Enrichment Prediction

### Abstract

By analyzing 34,524 probes in 10 cancer types across 20 data sets (tumor versus
normal) from the Affymetrix GeneChip Human Genome U133 Plus 2.0 Array (GPL570)
microarray platform, we built a multioutput logistic regression classifier to
predict gene ontology (GO) enrichment for a small subset of important terms.
By using a leave-one-out cross validation classification accuracy, we displayed
that such a prediction is possible with accuracies ranging from as 45% - 80% for
22  gene ontology terms.  After implementing a forward feature selection on our
gene feature set accuracies improved, ranging from 80% - 100%. Our findings
potentially highlight that specific gene expression levels play a vital role in
mapping to gene ontology terms. This would also indicate that selected genes are
most responsible for the functional properties of a given gene product.

### Folder Structure

```
.
├── README.md                       [Project description and installation instructions]
├── config/config.yaml              [Configuration file for project]
├── docs
│   ├── images                      [Contains images used in the report and presentation]
│   ├── includes                    [Contains files used to generate the report]
│   ├── presentation.qmd            [Presentation in quarto format]
│   ├── presentation.pdf            [Presentation]
│   ├── report.pdf                  [Project report]
│   └── report.tex                  [Project report in LaTeX format]
├── index.yaml                      [Contains project details]
├── models                          [Contains trained models]
├── notebooks
│   ├── figures/                    [Contains figures generated by notebooks]
│   ├── project.ipynb               [Notebook used for analysis]
│   └── preliminaryanalysis.ipynb   [Notebook used for preliminary analysis]
└── scripts
    ├── cleandata.R                 [Prepares dataframes with GO terms and log fold changes]
    ├── loaddatasets.R              [Downloads datasets and performs DGE analysis]
    ├── loadmetadata.R              [Loads metadata for datasets]
    ├── prepdata.py                 [Prepares training, validation, and test data]
    └── process_batches.sh          [Runs all scripts to download and process data]
```

### Installation & Use

- All dependencies for this project are contained in `requirements.txt`. To
  install dependencies, run `pip3 install -r requirements.txt`
- To download the datasets, you can run `bash scripts/process_batches.sh`
  - This will download the datasets and perform the DGE analysis and GO
    enrichment. It will then prepare the data for use in the deep learning
    model.
  - Note that this will take a long time to run, and will require >20GB of
    storage space.
  - Alternatively, you can download the data from my Dropbox, [here](https://www.dropbox.com/s/if2x86765uc5l1k/data.tar.gz\?dl\=1) or by using the commands below.
    - This will create a folder called data in the PATH_TO_DATA_FOLDER
      directory.
    - Update the paths in the [config.yaml](config/config.yaml) file to point
      to the data folder.

```bash
# Download preprocessed data from Dropbox
wget "https://www.dropbox.com/s/if2x86765uc5l1k/data.tar.gz\?dl\=1" -O data.tar.gz
tar -xvf data.tar.gz -C PATH_TO_DATA_FOLDER
```

# Datasets
- CuMiDa
  - citation
  - 21 datasets on GPL570
    - other platforms?
- CuRiDa? - more training data

# DGE
- R
  - limma - linear models
  - Emprical bayes moderation
  - GO enrichment
    - FDR = 0.05
# Landmark genes

- based on D-GEX paper (from their repo)
  - need the citation from the list


# Packages
- limma
- Go.db
- hgu133plus2.db



Datasets:
- Load all DEE2 datasets that have corresponding GEO accessions
- Filter datasets with only 1 sample out
- Total of 11334, including over 300000 SRA runs

Processing (loaddatasets)
- 12 separate batches, each with ~945 datasets
- Download dataset
- Perform DGE analysis
- Perform Gene Ontology Enrichment
- Return matrix with log fold changes and GO terms

Processing (cleandata)
- Convert the downloaded data into a matrix that can be read by python

Processing (prepdata)
- Combine data from each batch into single dataframe
- Select only genes that are common across all data sets (How many?)
- One hot encode the GO terms (labels)
- Select GO terms that are present in

Samples: (7130,)
Features: (48978,)
X: (7130, 48978)
Y: (7130,)

Total GO terms: 11130
Y: (7130, 11130)

X_train: (4563, 48978)
Y_train: (4563, 11130)
X_val: (1141, 48978)
Y_val: (1141, 11130)
X_test: (1426, 48978)
Y_test: (1426, 11130)
