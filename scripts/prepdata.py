"""
prepdata.py
-----------

This script prepares the training, validation, and test data for the
model. It reads .feather files from the processed data directory and
saves them as .npy files in the data directory.
"""

from sklearn.model_selection import train_test_split
from feather import read_dataframe
from tqdm.auto import tqdm
import pandas as pd
import numpy as np
import yaml
import glob


def main():
    """Main function."""

    # Load list of data files
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)
    files = glob.glob(config["path"]["processed"] + "/*.feather")

    # Load data frames
    data = []
    with tqdm(total=len(files), desc="Reading data frames") as pbar:
        for file in files:
            data.append(read_dataframe(file))
            pbar.update(1)

    # Concatenate data frames
    data = pd.concat(data, axis=0).dropna(axis=1, how="any")

    # Spit data into features and labels
    samples = data["GSE"].to_numpy()
    features = data.drop(["GSE", "GO_terms"], axis=1).columns.to_numpy()
    X = data.drop(["GSE", "GO_terms"], axis=1).to_numpy()
    Y = data["GO_terms"]

    # Report shapes
    print("Samples:", samples.shape)
    print("Features:", features.shape)
    print("X:", X.shape)
    print("Y:", Y.shape)

    # One-Hot encode GO terms
    Y = Y.str.get_dummies(sep=";")
    print("Total GO terms:", Y.shape[1])

    # Sort Y by GO term frequency
    Y = Y.reindex(Y.sum(axis=0).sort_values(ascending=False).index, axis=1)
    labels = Y.columns.to_numpy()
    Y = Y.to_numpy()
    print("Y:", Y.shape)

    # Split data into training, validation, and test sets
    X_train, X_test, Y_train, Y_test = train_test_split(
        X, Y, test_size=config["split"]["test"], random_state=config["seed"]
    )
    X_train, X_val, Y_train, Y_val = train_test_split(
        X_train,
        Y_train,
        test_size=config["split"]["val"],
        random_state=config["seed"],
    )

    # Report shapes
    print("X_train:", X_train.shape)
    print("Y_train:", Y_train.shape)
    print("X_val:", X_val.shape)
    print("Y_val:", Y_val.shape)
    print("X_test:", X_test.shape)
    print("Y_test:", Y_test.shape)

    # Save data
    print("Saving data...")
    np.savez_compressed(
        config["path"]["final"] + "/labels.npz",
        labels=labels,
        samples=samples,
        features=features,
    )
    np.savez_compressed(
        config["path"]["final"] + "/train.npz", X=X_train, Y=Y_train
    )
    np.savez_compressed(config["path"]["final"] + "/val.npz", X=X_val, Y=Y_val)
    np.savez_compressed(
        config["path"]["final"] + "/test.npz", X=X_test, Y=Y_test
    )

    print("Done.")


if __name__ == "__main__":
    main()
