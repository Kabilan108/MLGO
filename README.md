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
├── README.md           [Contains project description and installation instructions]
├── index.yml           [Contains project details]
├── presentation.pptx   [Project presentation]
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
- In order to perform the analysis and generate results, please run the
  `project.ipynb` notebook. Note that running the notebook will likely take
  several hours depending on your system specs.
