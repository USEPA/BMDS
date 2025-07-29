__version__ = "25.2a1"  # see docs/development for versioning

import os

import matplotlib as mpl

from .batch import BatchResponse, BatchSession  # noqa: F401
from .constants import DistType as ContinuousDistType  # noqa: F401
from .constants import Models, PriorClass, PriorDistribution  # noqa: F401
from .datasets import (  # noqa: F401
    ContinuousDataset,
    ContinuousIndividualDataset,
    DichotomousDataset,
    NestedDichotomousDataset,
)
from .models.multi_tumor import Multitumor  # noqa: F401
from .session import Session  # noqa: F401
from .types.continuous import ContinuousRiskType  # noqa: F401
from .types.dichotomous import DichotomousRiskType  # noqa: F401
from .utils import citation  # noqa: F401

# pick the matplotlib Agg backend if unspecified in configuration
mpl.use(os.environ.get("MPLBACKEND", "Agg"))
