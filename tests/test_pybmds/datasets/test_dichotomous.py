from copy import deepcopy

import pytest
from pydantic import ValidationError

import pybmds
from pybmds.datasets.dichotomous import DichotomousDatasetSchema

dummy2 = [1, 2]
dummy3 = [1, 2, 3]
dummy4 = [1, 2, 3, 4]
dummy3_floats = [0.1, 0.2, 0.3]


class TestDichotomousDataset:
    def test_validation(self):
        # these should be valid
        pybmds.DichotomousDataset(doses=dummy3, ns=dummy3, incidences=dummy3)
        # some data adjustments result in non-integer based counts
        pybmds.DichotomousDataset(doses=dummy3_floats, ns=dummy3_floats, incidences=dummy3_floats)
        # these should raise errors
        with pytest.raises((IndexError, ValueError)):
            # insufficient number of dose groups
            pybmds.DichotomousDataset(doses=dummy2, ns=dummy2, incidences=dummy2)
        with pytest.raises((IndexError, ValueError)):
            # different sized lists
            pybmds.DichotomousDataset(doses=dummy4, ns=dummy3, incidences=dummy3)
        with pytest.raises((IndexError, ValueError)):
            # incidence > n
            pybmds.DichotomousDataset(doses=dummy3, ns=[3, 3, 3], incidences=[3, 3, 4])
        with pytest.raises((IndexError, ValueError)):
            # zero in ns data
            pybmds.DichotomousDataset(doses=dummy3, ns=[0, 2, 3], incidences=dummy3)

    def test_metadata(self):
        ds = pybmds.DichotomousDataset(
            doses=[0, 1.96, 5.69, 29.75], ns=[75, 49, 50, 49], incidences=[5, 1, 3, 14]
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

    def test_dfile_outputs(self):
        ds = pybmds.DichotomousDataset(doses=dummy3, ns=[5, 5, 5], incidences=[0, 1, 2])
        dfile = ds.as_dfile()
        expected = "Dose Incidence NEGATIVE_RESPONSE\n1.000000 0 5\n2.000000 1 4\n3.000000 2 3"
        assert dfile == expected

    def test_dose_drops(self):
        ddataset = pybmds.DichotomousDataset(
            doses=list(reversed([0, 1.96, 5.69, 29.75])),
            ns=list(reversed([75, 49, 50, 49])),
            incidences=list(reversed([5, 1, 3, 14])),
        )
        assert (
            ddataset.as_dfile()
            == "Dose Incidence NEGATIVE_RESPONSE\n0.000000 5 70\n1.960000 1 48\n5.690000 3 47\n29.750000 14 35"
        )
        ddataset.drop_dose()
        assert (
            ddataset.as_dfile()
            == "Dose Incidence NEGATIVE_RESPONSE\n0.000000 5 70\n1.960000 1 48\n5.690000 3 47"
        )
        with pytest.raises(ValueError):
            ddataset.drop_dose()

    def test_serialize(self):
        ds1 = pybmds.DichotomousDataset(
            doses=[1, 2, 3, 4], ns=[1, 2, 3, 4], incidences=[1, 2, 3, 4], id=123, name="test"
        )

        # make sure serialize looks correct
        serialized = ds1.serialize()
        assert serialized.model_dump(exclude={"plotting"}) == {
            "dtype": "D",
            "metadata": {
                "id": 123,
                "name": "test",
                "dose_units": "",
                "response_units": "",
                "dose_name": "",
                "response_name": "",
            },
            "doses": [1.0, 2.0, 3.0, 4.0],
            "ns": [1, 2, 3, 4],
            "incidences": [1, 2, 3, 4],
        }

        # make sure we get the correct class back
        ds2 = serialized.deserialize()
        assert isinstance(ds2, pybmds.DichotomousDataset)

        # make sure we get the same result back after deserializing
        assert ds1.serialize().model_dump() == ds2.serialize().model_dump()


class TestDichotomousDatasetSchema:
    def test_schema(self, ddataset):
        # check that cycling through serialization returns the same
        v1 = ddataset.serialize().model_dump()
        v2 = DichotomousDatasetSchema.model_validate(v1).deserialize().serialize().model_dump()
        assert v1 == v2

        data = deepcopy(v1)
        data["ns"] = [1, 2]
        with pytest.raises(ValidationError, match="Length"):
            DichotomousDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data["doses"] = data["doses"][:-1]
        with pytest.raises(
            ValidationError,
            match="Length of doses, ns, and incidences are not the same",
        ):
            DichotomousDatasetSchema.model_validate(data)

        # check incidence > n
        data = deepcopy(v1)
        data.update(incidences=[75, 49, 50, 49], ns=[5, 1, 3, 14])
        with pytest.raises(ValidationError, match="Incidence > N"):
            DichotomousDatasetSchema.model_validate(data)

        # min groups
        data = deepcopy(v1)
        data.update(doses=[1, 2], ns=[1, 2], incidences=[1, 2])
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            DichotomousDatasetSchema.model_validate(data)
