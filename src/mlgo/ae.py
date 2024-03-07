import os
from typing import List, Optional, Tuple

import plotly.graph_objects as go

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
import torch.nn.functional as F
from keras.datasets import fashion_mnist
from PIL import Image
from sklearn.decomposition import PCA
from torch import Tensor, nn, optim
from torch.cuda.amp import GradScaler
from torch.optim.optimizer import Optimizer
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
from torchsummary import summary
from tqdm.auto import tqdm

from glob import glob
import re


# Checking if CUDA is available
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {DEVICE} device")

# define path to store latent space embeddings
LSPATH = "data/latent_space"
if not os.path.exists(LSPATH):
    os.makedirs(LSPATH)


def plot_components(
    components: np.ndarray, labels: np.ndarray, variance: Optional[np.ndarray] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot the components of the PCA or other dimensionality reduction algorithm.

    Parameters
    ----------
    components : np.ndarray
        The components of the dimensionality reduction algorithm.
    labels : np.ndarray
        The labels of the data.
    variance : np.ndarray, optional
        The explained variance of each component, by default None
    """

    fig, axis = plt.subplots(figsize=(4, 4))

    for label in set(labels):
        axis.scatter(
            components[labels == label, 0],
            components[labels == label, 1],
            label=label,
            alpha=0.2,
            s=2,
        )

    axis.legend(bbox_to_anchor=(1.1, 1.05))

    if variance is not None:
        axis.set_xlabel(f"PCA 1 ({variance[0]:.2f}%)")
        axis.set_ylabel(f"PCA 2 ({variance[1]:.2f}%)")

    return fig, axis


def collate_fn(batch: List[Tuple[Tensor, int]]) -> Tuple[Tensor, Tensor]:
    """
    Collate function for FashionMNIST dataloader.

    Combines images and labels into separate tensors
    """

    process = transforms.Compose(
        [transforms.ToTensor(), transforms.Pad(2)]  # add 2px padding to each side
    )

    # process images in batch
    images = torch.stack([process(img) for img, _ in batch])

    # get labels
    labels = torch.LongTensor([label for _, label in batch])

    return images, labels


def plot_latent_embedding(
    embedding: np.ndarray, labels: np.ndarray, title: str = "Latent Space"
) -> go.Figure:
    """
    Plot the latent space of the model.

    Parameters
    ----------
    embedding : np.ndarray
        The latent space embedding.
    labels : np.ndarray
        The labels of the data.
    title : str, optional
        The title of the plot, by default "Latent Space"
    """

    return None


def train(
    model: nn.Module,
    train_loader: DataLoader,
    optimizer: Optimizer,
    criterion: nn.Module,
    scaler: GradScaler,
    save_embedding: bool = False,
) -> float:
    """
    Run a single training epoch.

    Parameters
    ----------
    model : nn.Module
        PyTorch model to train
    train_loader : DataLoader
        data loader for training data
    optimizer : Optimizer
        training optimizer
    criterion : nn.Module
        loss function
    scaler : GradScaler
        scaler for mixed precision training
    save_embedding : bool, optional
        save latent space embedding, by default False

    Returns
    -------
    float
        average loss
    """

    # train model
    model.train()

    # keep track of losses
    train_loss = 0.0

    # tqdm bar for training batches
    train_bar = tqdm(
        total=len(train_loader),
        dynamic_ncols=True,
        leave=False,
        position=0,
        desc="Training",
    )

    # train model
    for batch, (images, labels) in enumerate(train_loader):
        # move to device
        images = images.to(DEVICE)

        # use mixed precision training
        with torch.cuda.amp.autocast():
            embedding = model.encoder(images)
            reconstruction = model.decoder(embedding)
            loss = criterion(reconstruction, images)

        # backprop
        optimizer.zero_grad()
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()

        # update loss
        train_loss += loss.item() * images.size(0)

        # update tqdm bar
        train_bar.set_postfix(
            loss=f"{train_loss / (batch + 1):.4f}",
            lr=f"{optimizer.param_groups[0]['lr']:.2e}",
        )
        train_bar.update()

        # save latent space embedding
        if save_embedding and batch % 10 == 0:
            # current epoch (add 1 because we start at 0)
            epoch = model.TRAINED_EPOCHS + 1

            # convert embedding and labels to numpy arrays
            embedding = embedding.detach().cpu().numpy()
            labels = labels.detach().cpu().numpy()

            # save embedding and labels
            np.savez(
                f"{LSPATH}/latent_space_E{epoch}_B{batch}.npz",
                points=embedding,
                labels=labels,
            )

    # update trained epochs
    model.TRAINED_EPOCHS += 1

    # close tqdm bar
    train_bar.close()

    # return average loss
    return train_loss / len(train_loader.dataset)


def validate(model: nn.Module, valid_loader: DataLoader, criterion: nn.Module) -> float:
    """
    Run a single validation epoch.

    Parameters
    ----------
    model : nn.Module
        PyTorch model to train
    valid_loader : DataLoader
        data loader for validation data
    criterion : nn.Module
        loss function

    Returns
    -------
    float
        average loss
    """

    # eval model
    model.eval()

    # keep track of losses
    valid_loss = 0.0

    # tqdm bar for training batches
    train_bar = tqdm(
        total=len(valid_loader),
        dynamic_ncols=True,
        leave=False,
        position=0,
        desc="Validation",
    )

    # train model
    for batch, (images, _) in enumerate(valid_loader):
        # move to device
        images = images.to(DEVICE)

        # forward pass
        with torch.no_grad():
            y_hat = model(images)
            loss = criterion(y_hat, images)

        # update loss
        valid_loss += loss.item() * images.size(0)

        # update tqdm bar
        train_bar.set_postfix(
            loss=f"{valid_loss / (batch + 1):.4f}",
        )
        train_bar.update()

    # close tqdm bar
    train_bar.close()

    # return average loss
    return valid_loss / len(valid_loader.dataset)
