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
├── README.md                       [Rroject description and installation instructions]
├── config                          [Configuration files]
├── docs
│   ├── _quarto.yml
│   ├── images
│   │   ├── cancer-metastasis.png
│   │   └── cumida.png
│   └── presentation.qmd
├── index.yaml
├── models
├── notebooks
│   ├── Model-Development.ipynb
│   ├── figures
│   └── preliminaryanalysis.ipynb
└── scripts
    ├── cleandata.R
    ├── loaddatasets.R
    ├── loadmetadata.R
    ├── prepdata.py
    └── process_batches.sh
```

```
.
├── README.md           [Contains project description and installation instructions]
├── index.yaml          [Contains project details]
├── presentation.pdf   [Project presentation]
├── report.docx         [Project report]
├── requirements.txt    [Contains list of project dependencies]
├── results             [Contains tables and figures generated]
│   ├── 2022065-logreg-before-feature-selection.png
│   ├── selected_datasets.tsv
│   ├── table1-cancer-counts.tsv
│   └── table2-selected-genes.tsv
├── src                 [Contains project source code]
│   ├── dataprep.py     [Script for preparing the machine learning dataset]
│   ├── project.ipynb   [Notebook for running analysis]
│   └── tools.py        [Module with custom functions and classes]
└── thumb.png           [Project thumbnail]
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
