import pytest

import pybmds


# dataset tests
@pytest.mark.mpl_image_compare
def test_cdataset_plot(cdataset):
    return cdataset.plot()


@pytest.mark.mpl_image_compare
def test_cidataset_plot(cidataset):
    return cidataset.plot()


@pytest.mark.mpl_image_compare
def test_ddataset_plot(ddataset):
    return ddataset.plot()


# test custom axes
@pytest.mark.mpl_image_compare
def test_cdataset_plot_customized():
    return pybmds.ContinuousDataset(
        doses=[0, 10, 50, 150, 400],
        ns=[111, 142, 143, 93, 42],
        means=[2.112, 2.095, 1.956, 1.587, 1.254],
        stdevs=[0.235, 0.209, 0.231, 0.263, 0.159],
        name="Example",
        dose_units="μg/m³",
        response_units="%",
        dose_name="Dose",
        response_name="Relative liver weight",
    ).plot()


@pytest.mark.mpl_image_compare
def test_cidataset_plot_customized():
    # fmt: off
    return pybmds.ContinuousIndividualDataset(
        doses=[
            0, 0, 0, 0, 0, 0, 0, 0,
            0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
            1, 1, 1, 1, 1, 1,
            10, 10, 10, 10, 10, 10,
            100, 100, 100, 100, 100, 100,
            300, 300, 300, 300, 300, 300,
            500, 500, 500, 500, 500, 500,
        ],
        responses=[
            8.1079, 9.3063, 9.7431, 9.7814, 10.0517, 10.6132, 10.7509, 11.0567,
            9.1556, 9.6821, 9.8256, 10.2095, 10.2222, 12.0382,
            9.5661, 9.7059, 9.9905, 10.2716, 10.471, 11.0602,
            8.8514, 10.0107, 10.0854, 10.5683, 11.1394, 11.4875,
            9.5427, 9.7211, 9.8267, 10.0231, 10.1833, 10.8685,
            10.368, 10.5176, 11.3168, 12.002, 12.1186, 12.6368,
            9.9572, 10.1347, 10.7743, 11.0571, 11.1564, 12.0368
        ],
        name="Example",
        dose_units="μg/m³",
        response_units="%",
        dose_name="Dose",
        response_name="Relative liver weight",
    ).plot()
    # fmt: on


@pytest.mark.mpl_image_compare
def test_ddataset_plot_customized():
    return pybmds.DichotomousDataset(
        doses=[0, 1.96, 5.69, 29.75],
        ns=[75, 49, 50, 49],
        incidences=[5, 1, 3, 14],
        name="Example",
        dose_units="μg/m³",
        response_units="%",
        dose_name="Dose",
        response_name="Fraction observed",
    ).plot()
