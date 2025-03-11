from copy import deepcopy

import pytest
from pydantic import ValidationError

import pybmds
from pybmds.datasets.continuous import ContinuousDatasetSchema, ContinuousIndividualDatasetSchema
from pybmds.datasets.dichotomous import DichotomousDatasetSchema


class TestBaseDatasetFunctionality:
    # we use dichotomous as subclass to test; these methods are on base

    def test_dose_space(self):
        ds = pybmds.DichotomousDataset(
            doses=[0, 1.96, 5.69, 29.75], ns=[75, 49, 50, 49], incidences=[5, 1, 3, 14]
        )
        xs = ds.dose_linspace()
        assert xs.min() == pybmds.constants.ZEROISH
        assert xs.max() == 29.75
        assert xs.size == 100

        # add extra values
        xs = ds.dose_linspace(extra_values=[25, 50], n=10)
        assert xs.min() == pybmds.constants.ZEROISH
        assert xs.max() == 50
        assert xs.size == 10

    def test_reporting_metadata(self):
        ds = pybmds.DichotomousDataset(
            doses=[0, 1.96, 5.69, 29.75], ns=[75, 49, 50, 49], incidences=[5, 1, 3, 14]
        )
        assert ds._get_dose_units_text() == ""
        assert ds._get_response_units_text() == ""
        assert ds._get_dataset_name() == "BMDS output results"

        ds = pybmds.DichotomousDataset(
            doses=[0, 1.96, 5.69, 29.75],
            ns=[75, 49, 50, 49],
            incidences=[5, 1, 3, 14],
            id=123,
            name="example dataset",
            dose_units="mg/kg/d",
            response_units="ug/m3",
            dose_name="Intravenous",
            response_name="Volume",
        )
        assert ds._get_dose_units_text() == " (mg/kg/d)"
        assert ds._get_response_units_text() == " (ug/m3)"
        assert ds._get_dataset_name() == "example dataset"

        ds = pybmds.DichotomousDataset(
            doses=[0, 1.96, 5.69, 29.75],
            ns=[75, 49, 50, 49],
            incidences=[5, 1, 3, 14],
            id=123,
        )
        assert ds._get_dataset_name() == "Dataset #123"

    def test_extra_metadata(self):
        # extra metadata is allowed, but not used in bmds
        ds = pybmds.DichotomousDataset(
            doses=[0, 1.96, 5.69, 29.75],
            ns=[75, 49, 50, 49],
            incidences=[5, 1, 3, 14],
            id=123,
            extra=[1, 2, 3],
        )
        assert ds.metadata.extra == [1, 2, 3]
        assert ds.metadata.model_dump()["extra"] == [1, 2, 3]


class TestSchema:
    def test_dichotomous(self, ddataset):
        # check that cycling through serialization returns the same
        v1 = ddataset.serialize().model_dump()
        v2 = DichotomousDatasetSchema.model_validate(v1).deserialize().serialize().model_dump()
        assert v1 == v2

        data = deepcopy(v1)
        data["ns"] = [1, 2]
        with pytest.raises(ValidationError, match="Length"):
            DichotomousDatasetSchema.model_validate(data)

        data = {
            "dtype": "D",
            "doses": [0, 10],
            "ns": [20, 20],
            "incidences": [0, 0],
            "metadata": {},
        }
        # but not ok for a standard dataset
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            DichotomousDatasetSchema.model_validate(data)

    def test_continuous(self, cdataset):
        # check that cycling through serialization returns the same
        v1 = cdataset.serialize().model_dump()
        v2 = ContinuousDatasetSchema.model_validate(v1).deserialize().serialize().model_dump()
        assert v1 == v2

        data = deepcopy(v1)
        data["ns"] = [1, 2]
        with pytest.raises(ValidationError, match="Length"):
            ContinuousDatasetSchema.model_validate(data)

        data = {
            "dtype": "C",
            "doses": [0, 10],
            "ns": [20, 20],
            "means": [0, 0],
            "stdevs": [1, 1],
            "metadata": {},
        }
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            ContinuousDatasetSchema.model_validate(data)

    def test_ci_continuous(self, cidataset):
        # check that cycling through serialization returns the same
        v1 = cidataset.serialize().model_dump()
        v2 = (
            ContinuousIndividualDatasetSchema.model_validate(v1)
            .deserialize()
            .serialize()
            .model_dump()
        )
        assert v1 == v2

        data = deepcopy(v1)
        data["doses"] = [1, 2]
        with pytest.raises(ValidationError, match="Length of doses and responses are not the same"):
            ContinuousIndividualDatasetSchema.model_validate(data)

        data = {
            "dtype": "CI",
            "doses": [0, 10],
            "responses": [20, 20],
            "metadata": {},
        }
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            ContinuousIndividualDatasetSchema.model_validate(data)
