import abc
from typing import TypeVar

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pydantic import BaseModel, ConfigDict

from .. import plotting
from ..constants import ZEROISH, Dtype


class DatasetMetadata(BaseModel):
    id: int | None = None
    name: str = ""
    dose_units: str = ""
    response_units: str = ""
    dose_name: str = ""
    response_name: str = ""
    model_config = ConfigDict(extra="allow")

    def get_name(self):
        if self.name:
            return self.name
        if self.id:
            return f"Dataset #{self.id}"
        return "BMDS output results"


class DatasetBase(abc.ABC):
    # Abstract parent-class for dataset-types.

    dtype: Dtype
    metadata: DatasetMetadata

    DEFAULT_XLABEL = "Dose"
    DEFAULT_YLABEL = "Response"

    @abc.abstractmethod
    def _validate(self): ...

    @abc.abstractmethod
    def as_dfile(self): ...

    @abc.abstractmethod
    def plot(self) -> Figure: ...

    def setup_plot(self, figsize: tuple[float, float] | None = None) -> Axes:
        """
        Return a matplotlib Axes of the dose-response dataset.
        """
        fig = plotting.create_empty_figure(figsize=figsize)
        ax = fig.gca()
        ax.set_xlabel(self.get_xlabel())
        ax.set_ylabel(self.get_ylabel())
        ax.margins(plotting.PLOT_MARGINS)
        ax.set_title(self._get_dataset_name())

        # set x bounds based on input data
        min_x = np.min(self.doses)
        max_x = np.max(self.doses)
        x_range = max_x - min_x
        ax.set_xlim(min_x - 0.05 * x_range, max_x + 0.05 * x_range)

        return ax

    @abc.abstractmethod
    def drop_dose(self): ...

    @property
    def num_dose_groups(self):
        return len(set(self.doses))

    def to_dict(self):
        return self.serialize().model_dump()

    @property
    def dose_linspace(self) -> np.ndarray:
        if not hasattr(self, "_dose_linspace"):
            self._dose_linspace = np.linspace(np.min(self.doses), np.max(self.doses), 100)
            self._dose_linspace[self._dose_linspace == 0] = ZEROISH
        return self._dose_linspace

    def _get_dose_units_text(self) -> str:
        if self.metadata.dose_units:
            return f" ({self.metadata.dose_units})"
        return ""

    def _get_response_units_text(self) -> str:
        if self.metadata.response_units:
            return f" ({self.metadata.response_units})"
        return ""

    def _get_dataset_name(self) -> str:
        return self.metadata.get_name()

    def get_xlabel(self):
        label = self.DEFAULT_XLABEL
        if self.metadata.dose_name:
            label = self.metadata.dose_name
        if self.metadata.dose_units:
            label += f" ({self.metadata.dose_units})"
        return label

    def get_ylabel(self):
        label = self.DEFAULT_YLABEL
        if self.metadata.response_name:
            label = self.metadata.response_name
        if self.metadata.response_units:
            label += f" ({self.metadata.response_units})"
        return label

    @abc.abstractmethod
    def serialize(self) -> "DatasetSchemaBase": ...

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            dataset_id=self.metadata.id,
            dataset_name=self.metadata.name,
            dataset_dose_name=self.metadata.dose_name,
            dataset_dose_units=self.metadata.dose_units,
            dataset_response_name=self.metadata.response_name,
            dataset_response_units=self.metadata.response_units,
        )

    @abc.abstractmethod
    def rows(self, extras: dict) -> list[dict]:
        """Return a list of rows; one for each item in a dataset"""
        ...


DatasetType = TypeVar("DatasetType", bound=DatasetBase)


class DatasetSchemaBase(BaseModel, abc.ABC):
    @classmethod
    def get_subclass(cls, dtype: Dtype) -> BaseModel:
        from .continuous import ContinuousDatasetSchema, ContinuousIndividualDatasetSchema
        from .dichotomous import DichotomousDatasetSchema
        from .nested_dichotomous import NestedDichotomousDatasetSchema

        _dataset_schema_map: dict = {
            Dtype.CONTINUOUS: ContinuousDatasetSchema,
            Dtype.CONTINUOUS_INDIVIDUAL: ContinuousIndividualDatasetSchema,
            Dtype.DICHOTOMOUS: DichotomousDatasetSchema,
            Dtype.NESTED_DICHOTOMOUS: NestedDichotomousDatasetSchema,
        }
        try:
            return _dataset_schema_map[dtype]
        except KeyError:
            raise ValueError(f"Unknown dtype: {dtype}") from None

    @abc.abstractmethod
    def deserialize(self) -> DatasetBase: ...


class DatasetPlottingSchema(BaseModel):
    mean: list[float] | None = None
    ll: list[float]
    ul: list[float]
