"""
utils.py
--------

Utility functions for the MLGO project.
"""

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from numpy import ndarray

import matplotlib.pyplot as plt
import seaborn as sns

from typing import List, Tuple


def get_layer_sizes(num_layers: int, total_units: int) -> List[int]:
    """
    Get the number of units in each layer of a neural network based on the
    total number of units and the number of layers. Ensures that subsequent
    layers have roughly half the number of units as the previous layer.

    This is based on the series sum:
    scale * (1 + 1/2 + 1/4 + ... + 1/2^(num_layers - 1)) = scale * (2 - 1/2^num_layers)
                                                         = total_units

    Parameters
    ----------
    num_layers : int
        Number of layers in the network
    total_units : int
        Total number of units in the network

    Returns
    -------
    List[int]
        List of the number of units in each layer
    """

    layer_sizes = []
    scale = total_units / (2 - 0.5 ** (num_layers - 1))
    for i in range(num_layers):
        layer_sizes.append(int(round(scale / 2**i)))
    return layer_sizes


def build_autoencoder(num_layers: int, total_units: int, input_dim: int) -> Sequential:
    """
    Build an autoencoder network.

    This will return an autoencoder neural network with the specified number
    of layers and total number of units. The number of units in each layer
    will be determined by the `get_layer_sizes` function.

    Parameters
    ----------
    num_layers : int
        Number of layers in the network
    total_units : int
        Total number of units in the network
    input_dim : int
        Number of input features

    Returns
    -------
    Sequential
        Keras Sequential model
    """

    # Get layer sizes
    layer_sizes = get_layer_sizes(num_layers, total_units)

    # Define the encoder network
    encoder = Sequential()
    for i, units in enumerate(layer_sizes):
        if i == 0:
            encoder.add(Dense(units, activation="relu", input_shape=(input_dim,)))
        else:
            encoder.add(Dense(units, activation="relu"))

    # Define the decoder network
    decoder = Sequential()
    for units in reversed(layer_sizes[:-1]):
        decoder.add(Dense(units, activation="relu"))
    decoder.add(Dense(input_dim, activation="sigmoid"))

    # Combine the encoder and decoder networks
    autoencoder = Sequential([encoder, decoder])

    return autoencoder


def plotpca(Y: ndarray, V: ndarray, samples: ndarray = None) -> Tuple[Figure, Axes]:
    """
    Create a 2D PCA plot.

    Parameters
    ----------
    Y : ndarray
        PCA data with >= 2 dimensions
    V : ndarray
        Explained variance ratio
    samples : ndarray
        Sample labels

    Returns
    -------
    Tuple[Figure, Axes]
        Matplotlib Figure and Axes objects
    """

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    g = sns.scatterplot(
        x=Y[:, 0],
        y=Y[:, 1],
        hue=samples,
        s=150,
        alpha=0.5,
        legend="full",
        palette="Dark2",
        ax=ax,
    )

    g.legend(loc="upper right", fontsize=12, frameon=False)

    ax.set_xlabel(f"PC1 ({V[0]:.2%})", fontsize=14)
    ax.set_ylabel(f"PC2 ({V[1]:.2%})", fontsize=14)
    ax.set_title("PCA", fontsize=16)

    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)
    ax.tick_params(axis="both", labelsize=12)
    ax.minorticks_on()
    ax.grid(which="minor", linestyle=":", linewidth="0.5", color="black", alpha=0.4)
    ax.grid(which="major", linestyle="-", linewidth="0.5", color="black", alpha=0.7)
    ax.axis("square")

    return fig, ax
