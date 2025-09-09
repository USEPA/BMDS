__version__ = "25.2a2"  # see docs/development for versioning

import os

import matplotlib as mpl

from .batch import BatchResponse, BatchSession
from .constants import DistType as ContinuousDistType
from .constants import Models, PriorClass, PriorDistribution
from .datasets import (
    ContinuousDataset,
    ContinuousIndividualDataset,
    DichotomousDataset,
    NestedDichotomousDataset,
)
from .models.multi_tumor import Multitumor
from .session import Session
from .types.continuous import ContinuousRiskType
from .types.dichotomous import DichotomousRiskType
from .utils import citation

__all__ = [
    "BatchResponse",
    "BatchSession",
    "ContinuousDataset",
    "ContinuousDistType",
    "ContinuousIndividualDataset",
    "ContinuousRiskType",
    "DichotomousDataset",
    "DichotomousRiskType",
    "Models",
    "Multitumor",
    "NestedDichotomousDataset",
    "PriorClass",
    "PriorDistribution",
    "Session",
    "citation",
]

# pick the matplotlib Agg backend if unspecified in configuration
mpl.use(os.environ.get("MPLBACKEND", "Agg"))
