from copy import deepcopy

import numpy as np
import pytest
from pydantic import ValidationError

import pybmds
from pybmds.datasets.continuous import ContinuousDatasetSchema, ContinuousIndividualDatasetSchema

dummy2 = [1, 2]
dummy3 = [1, 2, 3]
dummy4 = [1, 2, 3, 4]
dummy3_dups = [0, 0, 1]
dummy3_floats = [0.1, 0.2, 0.3]


class TestContinuousSummaryDataset:
    def test_validation(self):
        # these should be valid
        pybmds.ContinuousDataset(doses=dummy3, ns=dummy3, means=dummy3, stdevs=dummy3)
        pybmds.ContinuousDataset(doses=dummy3, ns=dummy3_floats, means=dummy3, stdevs=dummy3)

        # these should raise exceptions
        with pytest.raises(ValueError, match="Must have 3 or more dose groups"):
            pybmds.ContinuousDataset(doses=dummy2, ns=dummy2, means=dummy2, stdevs=dummy2)
        with pytest.raises(IndexError):
            pybmds.ContinuousDataset(doses=dummy4, ns=dummy3, means=dummy3, stdevs=dummy3)
        with pytest.raises(ValueError, match="N must be positive and non-zero"):
            pybmds.ContinuousDataset(doses=dummy3, ns=[0, 2, 3], means=dummy3, stdevs=dummy3)
        with pytest.raises(ValueError, match="Stdev must be greater than zero"):
            pybmds.ContinuousDataset(doses=dummy3, ns=dummy3, means=dummy3, stdevs=[0, 1, 2])

    def test_extra_kwargs(self):
        ds = pybmds.ContinuousDataset(
            doses=[0, 10, 50, 150, 400],
            ns=[111, 142, 143, 93, 42],
            means=[2.112, 2.095, 1.956, 1.587, 1.254],
            stdevs=[0.235, 0.209, 0.231, 0.263, 0.159],
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
        assert ds.get_ylabel() == "Response"

        ds = pybmds.ContinuousDataset(
            doses=[0, 10, 50, 150, 400],
            ns=[111, 142, 143, 93, 42],
            means=[2.112, 2.095, 1.956, 1.587, 1.254],
            stdevs=[0.235, 0.209, 0.231, 0.263, 0.159],
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

    def test_is_increasing(self):
        ds = pybmds.ContinuousDataset(doses=dummy4, ns=dummy4, means=dummy4, stdevs=dummy4)
        assert ds.is_increasing is True

        rev = list(reversed(dummy4))
        ds = pybmds.ContinuousDataset(doses=dummy4, ns=dummy4, means=rev, stdevs=dummy4)
        assert ds.is_increasing is False

        ds = pybmds.ContinuousDataset(doses=dummy4, ns=dummy4, means=[1, 2, 3, 0], stdevs=dummy4)
        assert ds.is_increasing is True

        ds = pybmds.ContinuousDataset(doses=dummy4, ns=dummy4, means=[1, 3, 2, 1], stdevs=dummy4)
        assert ds.is_increasing is True

        ds = pybmds.ContinuousDataset(doses=dummy4, ns=dummy4, means=[0, 2, -1, 0], stdevs=dummy4)
        assert ds.is_increasing is True

    def test_anova(self, cdataset):
        actual = cdataset.get_anova_report()
        expected = "                     Tests of Interest    \n   Test    -2*log(Likelihood Ratio)  Test df        P-Value    \n   Test 1              475.068          8                0\n   Test 2              15.4005          4         0.003939\n   Test 3              15.4005          4         0.003939"
        assert actual == expected

    def test_dfile_outputs(self):
        ds = pybmds.ContinuousDataset(doses=dummy3, ns=dummy3, means=dummy3, stdevs=dummy3)
        dfile = ds.as_dfile()
        expected = "Dose NumAnimals Response Stdev\n1.000000 1 1.000000 1.000000\n2.000000 2 2.000000 2.000000\n3.000000 3 3.000000 3.000000"
        assert dfile == expected

    def test_dataset_reporting_options(self, cdataset):
        # test defaults
        assert cdataset._get_dose_units_text() == ""
        assert cdataset._get_response_units_text() == ""
        assert cdataset._get_dataset_name() == "BMDS output results"

        # test settings
        ds = pybmds.ContinuousDataset(
            doses=[0, 10, 50, 150, 400],
            ns=[111, 142, 143, 93, 42],
            means=[2.112, 2.095, 1.956, 1.587, 1.254],
            stdevs=[0.235, 0.209, 0.231, 0.263, 0.159],
            dose_units="μg/m³",
            response_units="mg/kg",
            name="Smith 2017: Relative Liver Weight in Male SD Rats",
        )
        assert ds._get_dataset_name() == "Smith 2017: Relative Liver Weight in Male SD Rats"
        assert ds._get_dose_units_text() == " (μg/m³)"
        assert ds._get_response_units_text() == " (mg/kg)"

    def test_dose_drops(self, cidataset):
        cdataset = pybmds.ContinuousDataset(
            doses=list(reversed([0, 10, 50, 150, 400])),
            ns=list(reversed([111, 142, 143, 93, 42])),
            means=list(reversed([2.112, 2.095, 1.956, 1.587, 1.254])),
            stdevs=list(reversed([0.235, 0.209, 0.231, 0.263, 0.159])),
        )

        assert (
            cdataset.as_dfile()
            == "Dose NumAnimals Response Stdev\n0.000000 111 2.112000 0.235000\n10.000000 142 2.095000 0.209000\n50.000000 143 1.956000 0.231000\n150.000000 93 1.587000 0.263000\n400.000000 42 1.254000 0.159000"
        )
        cdataset.drop_dose()
        assert (
            cdataset.as_dfile()
            == "Dose NumAnimals Response Stdev\n0.000000 111 2.112000 0.235000\n10.000000 142 2.095000 0.209000\n50.000000 143 1.956000 0.231000\n150.000000 93 1.587000 0.263000"
        )
        cdataset.drop_dose()
        assert (
            cdataset.as_dfile()
            == "Dose NumAnimals Response Stdev\n0.000000 111 2.112000 0.235000\n10.000000 142 2.095000 0.209000\n50.000000 143 1.956000 0.231000"
        )
        with pytest.raises(ValueError):
            cdataset.drop_dose()

    def test_serialize(self):
        ds1 = pybmds.ContinuousDataset(
            doses=[1, 2, 3, 4],
            ns=[1, 2, 3, 4],
            means=[1, 2, 3, 4],
            stdevs=[1, 2, 3, 4],
            id=123,
            name="test",
        )

        # make sure serialize looks correct
        serialized = ds1.serialize()
        assert serialized.model_dump(exclude={"plotting"}) == {
            "dtype": "C",
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
            "means": [1.0, 2.0, 3.0, 4.0],
            "stdevs": [1.0, 2.0, 3.0, 4.0],
            "anova": None,
        }

        # make sure we get the correct class back
        ds2 = serialized.deserialize()
        assert isinstance(ds2, pybmds.ContinuousDataset)
        assert not isinstance(ds2, pybmds.ContinuousIndividualDataset)

        # make sure we get the same result back after deserializing
        assert ds1.serialize().model_dump() == ds2.serialize().model_dump()

    def test_simulate_individual_dataset(self):
        # test a simulated dataset is created fits desired distribution
        dataset = pybmds.ContinuousDataset(
            doses=[0, 1, 5], ns=[10, 9, 8], means=[1, 2, 3], stdevs=[0.1, 0.2, 0.3]
        )
        fabricated = dataset.simulate_individual_dataset()
        assert isinstance(fabricated, pybmds.ContinuousIndividualDataset)
        assert len(fabricated.responses) == sum(dataset.ns)
        assert len(fabricated.individual_doses) == sum(dataset.ns)
        assert np.allclose(dataset.doses, fabricated.doses)
        assert np.allclose(dataset.ns, fabricated.ns)
        assert np.allclose(dataset.means, fabricated.means)
        assert np.allclose(dataset.stdevs, fabricated.stdevs)

    def test_simulate_individual_dataset_failure(self):
        # should fail with large SD close to zero w/ impose_positivity
        dataset = pybmds.ContinuousDataset(
            doses=[1, 2, 3], ns=[10] * 3, means=[1] * 3, stdevs=[5] * 3
        )
        with pytest.raises(ValueError, match="Could not generate a valid sample"):
            dataset.simulate_individual_dataset(
                seed=42, tolerance=1e-10, max_iterations=1, impose_positivity=True
            )


class TestContinuousIndividualDataset:
    def test_validation(self):
        # these should be valid
        pybmds.ContinuousIndividualDataset(
            doses=[1, 1, 2, 2, 3, 3, 4, 4], responses=[1, 2, 3, 4, 5, 6, 7, 8]
        )
        # different sized lists
        with pytest.raises(ValueError, match="All arrays must be of the same length"):
            pybmds.ContinuousIndividualDataset(doses=[1, 2, 3, 4], responses=[1, 2, 3])
        # < 3 dose groups
        with pytest.raises(ValueError, match="Must have 3 or more dose groups"):
            pybmds.ContinuousIndividualDataset(doses=[1, 2], responses=[1, 2])
        # n <= 1
        with pytest.raises(ValueError, match="more than one response for each dose group"):
            pybmds.ContinuousIndividualDataset(doses=[1, 2, 3], responses=[1, 2, 3])

    def test_extra_kwargs(self):
        ds = pybmds.ContinuousIndividualDataset(
            doses=[0, 0, 1, 1, 2, 2, 3, 3],
            responses=[8.1079, 9.3063, 9.7431, 9.7814, 10.0517, 10.6132, 10.7509, 11.0567],
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
        assert ds.get_ylabel() == "Response"

        ds = pybmds.ContinuousIndividualDataset(
            doses=[0, 0, 1, 1, 2, 2, 3, 3],
            responses=[8.1079, 9.3063, 9.7431, 9.7814, 10.0517, 10.6132, 10.7509, 11.0567],
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
        ds = pybmds.ContinuousIndividualDataset(
            doses=[1, 1, 2, 2, 3, 3], responses=[1, 2, 3, 4, 5, 6]
        )
        dfile = ds.as_dfile()
        expected = "Dose Response\n1.000000 1.000000\n1.000000 2.000000\n2.000000 3.000000\n2.000000 4.000000\n3.000000 5.000000\n3.000000 6.000000"
        assert dfile == expected

    def test_ci_summary_stats(self, cidataset):
        assert len(cidataset.doses) == 7
        assert np.allclose(cidataset.ns, [8, 6, 6, 6, 6, 6, 6])
        assert np.allclose(
            cidataset.means, [9.9264, 10.1889, 10.17755, 10.3571, 10.0275, 11.4933, 10.85275]
        )
        assert np.allclose(
            cidataset.stdevs, [0.87969, 0.90166, 0.50089, 0.85590, 0.42833, 0.83734, 0.690373]
        )

    def test_dose_drops(self, cidataset):
        assert (
            cidataset.as_dfile()
            == "Dose Response\n0.000000 8.107900\n0.000000 9.306300\n0.000000 9.743100\n0.000000 9.781400\n0.000000 10.051700\n0.000000 10.613200\n0.000000 10.750900\n0.000000 11.056700\n0.100000 9.155600\n0.100000 9.682100\n0.100000 9.825600\n0.100000 10.209500\n0.100000 10.222200\n0.100000 12.038200\n1.000000 9.566100\n1.000000 9.705900\n1.000000 9.990500\n1.000000 10.271600\n1.000000 10.471000\n1.000000 11.060200\n10.000000 8.851400\n10.000000 10.010700\n10.000000 10.085400\n10.000000 10.568300\n10.000000 11.139400\n10.000000 11.487500\n100.000000 9.542700\n100.000000 9.721100\n100.000000 9.826700\n100.000000 10.023100\n100.000000 10.183300\n100.000000 10.868500\n300.000000 10.368000\n300.000000 10.517600\n300.000000 11.316800\n300.000000 12.002000\n300.000000 12.118600\n300.000000 12.636800\n500.000000 9.957200\n500.000000 10.134700\n500.000000 10.774300\n500.000000 11.057100\n500.000000 11.156400\n500.000000 12.036800"
        )
        cidataset.drop_dose()
        assert (
            cidataset.as_dfile()
            == "Dose Response\n0.000000 8.107900\n0.000000 9.306300\n0.000000 9.743100\n0.000000 9.781400\n0.000000 10.051700\n0.000000 10.613200\n0.000000 10.750900\n0.000000 11.056700\n0.100000 9.155600\n0.100000 9.682100\n0.100000 9.825600\n0.100000 10.209500\n0.100000 10.222200\n0.100000 12.038200\n1.000000 9.566100\n1.000000 9.705900\n1.000000 9.990500\n1.000000 10.271600\n1.000000 10.471000\n1.000000 11.060200\n10.000000 8.851400\n10.000000 10.010700\n10.000000 10.085400\n10.000000 10.568300\n10.000000 11.139400\n10.000000 11.487500\n100.000000 9.542700\n100.000000 9.721100\n100.000000 9.826700\n100.000000 10.023100\n100.000000 10.183300\n100.000000 10.868500\n300.000000 10.368000\n300.000000 10.517600\n300.000000 11.316800\n300.000000 12.002000\n300.000000 12.118600\n300.000000 12.636800"
        )
        cidataset.drop_dose()
        assert (
            cidataset.as_dfile()
            == "Dose Response\n0.000000 8.107900\n0.000000 9.306300\n0.000000 9.743100\n0.000000 9.781400\n0.000000 10.051700\n0.000000 10.613200\n0.000000 10.750900\n0.000000 11.056700\n0.100000 9.155600\n0.100000 9.682100\n0.100000 9.825600\n0.100000 10.209500\n0.100000 10.222200\n0.100000 12.038200\n1.000000 9.566100\n1.000000 9.705900\n1.000000 9.990500\n1.000000 10.271600\n1.000000 10.471000\n1.000000 11.060200\n10.000000 8.851400\n10.000000 10.010700\n10.000000 10.085400\n10.000000 10.568300\n10.000000 11.139400\n10.000000 11.487500\n100.000000 9.542700\n100.000000 9.721100\n100.000000 9.826700\n100.000000 10.023100\n100.000000 10.183300\n100.000000 10.868500"
        )
        cidataset.drop_dose()
        cidataset.drop_dose()
        with pytest.raises(ValueError):
            cidataset.drop_dose()

    def test_serialize(self):
        ds1 = pybmds.ContinuousIndividualDataset(
            doses=[0, 0, 1, 1, 2, 2, 3, 3],
            responses=[8.1079, 9.3063, 9.7431, 9.7814, 10.0517, 10.6132, 10.7509, 11.0567],
            id=123,
            name="test",
        )

        # make sure serialize looks correct
        serialized = ds1.serialize()
        assert serialized.model_dump(exclude={"anova"}) == {
            "dtype": "CI",
            "metadata": {
                "id": 123,
                "name": "test",
                "dose_units": "",
                "response_units": "",
                "dose_name": "",
                "response_name": "",
            },
            "doses": [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
            "responses": [8.1079, 9.3063, 9.7431, 9.7814, 10.0517, 10.6132, 10.7509, 11.0567],
        }

        # make sure we get the correct class back
        ds2 = serialized.deserialize()
        assert not isinstance(ds2, pybmds.ContinuousDataset)
        assert isinstance(ds2, pybmds.ContinuousIndividualDataset)

        # make sure we get the same result back after deserializing
        assert ds1.serialize().model_dump() == ds2.serialize().model_dump()


class TestContinuousDatasetSchema:
    def test_schema(self, cdataset):
        # check that cycling through serialization returns the same
        v1 = cdataset.serialize().model_dump()
        v2 = ContinuousDatasetSchema.model_validate(v1).deserialize().serialize().model_dump()
        assert v1 == v2

        data = deepcopy(v1)
        data["doses"] = data["doses"][:-1]
        with pytest.raises(
            ValidationError,
            match="Length of doses, ns, means, and stdevs are not the same",
        ):
            ContinuousDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data.update(ns=[1, 2])
        with pytest.raises(ValidationError, match="Length"):
            ContinuousDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data.update(doses=[0, 10], ns=[20, 20], means=[0, 0], stdevs=[1, 1])
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            ContinuousDatasetSchema.model_validate(data)


class TestContinuousIndividualDatasetSchema:
    def test_schema(self, cidataset):
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
        data["doses"] = data["doses"][:-1]
        with pytest.raises(
            ValidationError,
            match="Length of doses and responses are not the same",
        ):
            ContinuousIndividualDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data["doses"] = [1, 2]
        with pytest.raises(ValidationError, match="Length of doses and responses are not the same"):
            ContinuousIndividualDatasetSchema.model_validate(data)

        data = deepcopy(v1)
        data.update(doses=[0, 10], responses=[20, 20])
        with pytest.raises(ValidationError, match="At least 3 groups are required"):
            ContinuousIndividualDatasetSchema.model_validate(data)
