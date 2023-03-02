"""
downloadgse.py

Download the 21 Gene Expression datasets from the Curated Microarray Database
and Prepare them for Gene Enrichment Analysis.
"""

# Imports from third party packages
from infoml.binf.data import CuMiDa
from infoml import CONFIG
from tqdm import tqdm
import yaml

# Imports from standard library
import shutil

# Set data and temp directories
CONFIG.datadir("data/raw")
CONFIG.tempdir("MLGO/cumida")


def main():
    """Main function"""

    # Load configuration
    with open("config/config.yaml") as f:
        conf = yaml.safe_load(f)["scripts"]["downloadgse.py"]

    # Download the data
    cumida = CuMiDa(CONFIG.tempdir())
    selected = cumida.index.query(conf["query"]).index.tolist()
    cumida.download(selected)

    # Prepare data for DGE analysis
    with tqdm(total=len(selected), desc="Preparing data") as pbar:
        for dataset in selected:
            fpath = CONFIG.datadir() / ("_".join(dataset) + ".tsv")

            gse = cumida.load(dataset).reset_index().drop(columns=["samples"])
            gse["type"] = gse["type"].apply(
                lambda x: "normal" if "normal" in x else "tumor"
            )
            gse = gse.set_index("type").T.reset_index(names=["Gene"])
            gse["Gene"] = gse["Gene"].str.split(".").str[0]

            gse.to_csv(fpath, sep="\t", index=False)

            pbar.update(1)

    # Remove temporary directory
    shutil.rmtree(CONFIG.tempdir())


if __name__ == "__main__":
    main()
