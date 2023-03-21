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
        datasets = []
        for dataset in selected:
            epath = CONFIG.datadir() / ("_".join(dataset) + "_exprs.tsv")
            dpath = CONFIG.datadir() / ("_".join(dataset) + "_data.tsv")

            gse = (
                cumida.load(dataset, probe_ids=True)
                .reset_index()
                .drop(columns=["samples"])
                .reset_index(names=["sample"])
            )
            gse["sample"] = "S" + (gse["sample"] + 1).astype(str)
            gse["type"] = gse["type"].apply(
                lambda x: "normal" if "normal" in x else "tumor"
            )

            exprs = (
                gse.drop(columns=["type"])
                .set_index("sample")
                .rename_axis("")
                .T
            )
            exprs.to_csv(epath, sep="\t", index=True)

            data = gse[["sample", "type"]].set_index("sample").rename_axis("")
            data.to_csv(dpath, sep="\t", index=True)

            datasets.append("_".join(dataset) + "\n")

            pbar.update(1)

    # Write list of datasets
    with open(CONFIG.datadir() / "datasets.txt", "w") as f:
        f.writelines(datasets)

    # Remove temporary directory
    shutil.rmtree(CONFIG.tempdir())


if __name__ == "__main__":
    main()
