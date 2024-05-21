from copy import deepcopy

import pytest
from pydantic import ValidationError

import pybmds
from pybmds.datasets.nested_dichotomous import NestedDichotomousDatasetSchema

dummy2 = [1, 2]
dummy3 = [1, 2, 3]
dummy4 = [1, 2, 3, 4]


class TestNestedDichotomousDataset:
    def test_validation(self):
        # these should be valid
        pybmds.NestedDichotomousDataset(
            doses=dummy3, litter_ns=dummy3, incidences=dummy3, litter_covariates=dummy3
        )
        # these should raise errors
        with pytest.raises((IndexError, ValueError)):
            # insufficient number of dose groups
            pybmds.NestedDichotomousDataset(
                doses=dummy2, litter_ns=dummy2, incidences=dummy2, litter_covariates=dummy2
            )
        with pytest.raises((IndexError, ValueError)):
            # different sized lists
            pybmds.NestedDichotomousDataset(
                doses=dummy4, litter_ns=dummy3, incidences=dummy3, litter_covariates=dummy3
            )

    def test_metadata(self):
        ds = pybmds.NestedDichotomousDataset(
            doses=dummy3, litter_ns=dummy3, incidences=dummy3, litter_covariates=dummy3
        )
        assert ds.to_dict()["metadata"] == {
            "id": None,
            "name": "",
            "dose_units": "",
            "response_units": "",
            "dose_name": "",
            "response_name": "",
        }

        assert ds.get_xlabel() == "Dose"
        assert ds.get_ylabel() == "Incidence"

        ds = pybmds.NestedDichotomousDataset(
            doses=dummy3,
            litter_ns=dummy3,
            incidences=dummy3,
            litter_covariates=dummy3,
            id=123,
            name="example dataset",
            dose_units="mg/kg/d",
            response_units="ug/m3",
            dose_name="Intravenous",
            response_name="Volume",
        )
        assert ds.to_dict()["metadata"] == {
            "id": 123,
            "name": "example dataset",
            "dose_units": "mg/kg/d",
            "response_units": "ug/m3",
            "dose_name": "Intravenous",
            "response_name": "Volume",
        }
        assert ds.get_xlabel() == "Intravenous (mg/kg/d)"
        assert ds.get_ylabel() == "Volume (ug/m3)"

    def test_drop_dose(self, nd_dataset):
        with pytest.raises(NotImplementedError):
            nd_dataset.drop_dose()

    def test_serialize(self, nd_dataset):
        ds1 = nd_dataset

        # make sure serialize looks correct
        # fmt: off
        assert ds1.serialize().model_dump() == {
            "dtype": "ND",
            "metadata": {
                "id": 123,
                "name": "",
                "dose_units": "",
                "response_units": "",
                "dose_name": "",
                "response_name": "",
            },
            "doses": [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                25, 25, 25, 25, 25, 25, 25, 25, 25,
                50, 50, 50, 50, 50, 50, 50, 50, 50,
            ],
            "litter_ns": [
                16, 9, 15, 14, 13, 9, 10, 14, 10, 11, 14,
                9, 14, 9, 13, 12, 10, 10, 11, 14,
                11, 11, 14, 11, 10, 11, 10, 15, 7,
            ],
            "incidences": [
                1, 1, 2, 3, 3, 0, 2, 2, 1, 2, 4,
                5, 6, 2, 6, 3, 1, 2, 4, 3,
                4, 5, 5, 4, 5, 4, 5, 6, 2,
            ],
            "litter_covariates": [
                16, 9, 15, 14, 13, 9, 10, 14, 10, 11, 14,
                9, 14, 9, 13, 12, 10, 10, 11, 14,
                11, 11, 14, 11, 10, 11, 10, 15, 7,
            ]

        }
        # fmt: on

        # make sure we get the correct class back
        ds2 = ds1.serialize().deserialize()
        assert isinstance(ds2, pybmds.NestedDichotomousDataset)

        # check serialization equality
        assert ds1.serialize().model_dump() == ds2.serialize().model_dump()


class TestNestedDichotomousSchema:
    def test_schema(self, nd_dataset):
        # check that cycling through serialization returns the same
        v1 = nd_dataset.serialize().model_dump()
        v2 = (
            NestedDichotomousDatasetSchema.model_validate(v1).deserialize().serialize().model_dump()
        )
        assert v1 == v2

        data = deepcopy(v1)
        data["doses"] = data["doses"][:-1]
        with pytest.raises(
            ValidationError,
            match="Length of dose, litter, incidence, and covariate are not the same",
        ):
            NestedDichotomousDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data.update(
            doses=[0, 10], litter_ns=[10, 10], incidences=[10, 10], litter_covariates=[1, 1]
        )
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            NestedDichotomousDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data["incidences"][0] = data["litter_ns"][0] + 1
        with pytest.raises(ValidationError, match="Incidence > N"):
            NestedDichotomousDatasetSchema.model_validate(data)
