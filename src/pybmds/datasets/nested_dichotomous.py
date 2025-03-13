import math
from typing import ClassVar

import numpy as np
from pydantic import model_validator

from .. import constants, plotting
from ..utils import str_list
from .base import DatasetBase, DatasetMetadata, DatasetSchemaBase


class NestedDichotomousDataset(DatasetBase):
    """
    Dataset object for nested dichotomous datasets.
    """

    _BMDS_DATASET_TYPE = 1  # group data
    MINIMUM_DOSE_GROUPS = 3
    dtype = constants.Dtype.NESTED_DICHOTOMOUS

    DEFAULT_YLABEL = "Incidence"

    def __init__(
        self,
        doses: list[float],
        litter_ns: list[float],
        incidences: list[float],
        litter_covariates: list[float],
        **metadata,
    ):
        self.doses = doses
        self.litter_ns = litter_ns
        self.incidences = incidences
        self.litter_covariates = litter_covariates
        self.metadata = DatasetMetadata.model_validate(metadata)
        self._sort_by_dose_group()
        self._validate()

    def _sort_by_dose_group(self):
        # use mergesort since it's a stable-sorting algorithm in numpy
        indexes = np.array(self.doses).argsort(kind="mergesort")
        for fld in ("doses", "litter_ns", "incidences", "litter_covariates"):
            arr = getattr(self, fld)
            setattr(self, fld, np.array(arr)[indexes].tolist())
        self._validate()

    def _validate(self):
        length = len(self.doses)
        if not all(
            len(lst) == length
            for lst in [self.doses, self.litter_ns, self.incidences, self.litter_covariates]
        ):
            raise ValueError("All input lists must be same length")

        if self.num_dose_groups < self.MINIMUM_DOSE_GROUPS:
            raise ValueError(
                f"Must have {self.MINIMUM_DOSE_GROUPS} or more dose groups after dropping doses"
            )

    def drop_dose(self):
        raise NotImplementedError("")

    def serialize(self) -> "NestedDichotomousDatasetSchema":
        return NestedDichotomousDatasetSchema(
            dtype=self.dtype,
            doses=self.doses,
            litter_ns=self.litter_ns,
            incidences=self.incidences,
            litter_covariates=self.litter_covariates,
            metadata=self.metadata,
        )

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        super().update_record(d)
        d.update(
            dataset_doses=str_list(self.doses),
            dataset_litter_ns=str_list(self.litter_ns),
            dataset_incidences=str_list(self.incidences),
            dataset_litter_covariates=str_list(self.litter_covariates),
        )

    def as_dfile(self):
        raise ValueError("N/A; requires BMDS3+ which doesn't use dfiles")

    def plot(self, figsize: tuple[float, float] | None = None):
        """
        Return a matplotlib figure of the dose-response dataset.

        Examples
        --------
        >>> fig = dataset.plot()
        >>> fig.show()
        >>> fig.clear()

        .. image:: ../tests/data/mpl/test_ddataset_plot.png
           :align: center
           :alt: Example generated BMD plot

        Returns
        -------
        out : matplotlib.figure.Figure
            A matplotlib figure representation of the dataset.
        """
        ax = self.setup_plot(figsize=figsize)
        ys = np.array(self.incidences) / np.array(self.litter_ns)
        ax.scatter(
            self.doses,
            ys,
            label="Fraction Affected",
            **plotting.DATASET_INDIVIDUAL_FORMAT,
        )
        ax.legend(**plotting.LEGEND_OPTS)
        return ax.get_figure()

    def rows(self, extras: dict | None = None) -> list[dict]:
        """Return a list of rows; one for each item in a dataset"""
        extra = self.metadata.model_dump()
        extra.update(extras or {})
        rows = []
        for dose, n, incidence, litter_covariate in zip(
            self.doses, self.litter_ns, self.incidences, self.litter_covariates, strict=True
        ):
            rows.append(
                {
                    **extra,
                    **dict(dose=dose, n=n, incidence=incidence, litter_covariate=litter_covariate),
                }
            )
        return rows

    def smean(self):
        return np.sum(self.litter_covariates) / float(len(self.litter_covariates))


class NestedDichotomousDatasetSchema(DatasetSchemaBase):
    dtype: constants.Dtype = constants.Dtype.NESTED_DICHOTOMOUS
    metadata: DatasetMetadata
    doses: list[float]
    litter_ns: list[int]
    incidences: list[int]
    litter_covariates: list[float]

    MIN_N: ClassVar = 3
    MAX_N: ClassVar = math.inf

    def deserialize(self) -> NestedDichotomousDataset:
        ds = NestedDichotomousDataset(
            doses=self.doses,
            litter_ns=self.litter_ns,
            incidences=self.incidences,
            litter_covariates=self.litter_covariates,
            **self.metadata.model_dump(),
        )
        return ds

    @model_validator(mode="after")
    def check_num_groups(self):
        n_doses = len(self.doses)
        n_litter_ns = len(self.litter_ns)
        n_incidences = len(self.incidences)
        n_litter_covariates = len(self.litter_covariates)
        if len(set([n_doses, n_litter_ns, n_incidences, n_litter_covariates])) > 1:
            raise ValueError("Length of dose, litter, incidence, and covariate are not the same")
        if n_doses < self.MIN_N:
            raise ValueError(f"At least {self.MIN_N} groups are required")
        if n_doses > self.MAX_N:
            raise ValueError(f"A maximum of {self.MAX_N} groups are allowed")
        for incidence, n in zip(self.incidences, self.litter_ns, strict=True):
            if incidence > n:
                raise ValueError(f"Incidence > N ({incidence} > {n})")
        return self
