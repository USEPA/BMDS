import json

import pytest

import pybmds
from pybmds.constants import BMDS_BLANK_VALUE, DistType, Models, PriorClass
from pybmds.exceptions import ConfigurationException
from pybmds.models import continuous
from pybmds.types.continuous import ContinuousModelSettings


class TestPriorOverrides:
    def test_poly(self, cdataset2, negative_cdataset):
        # fmt: off
        for settings, priors in [
            ({"disttype": DistType.normal, "is_increasing": True}, (0, 1e6)),
            ({"disttype": DistType.normal, "is_increasing": False}, (-1e6, 0)),
            ({"disttype": DistType.normal_ncv, "is_increasing": True}, (0, 18)),
            ({"disttype": DistType.normal_ncv, "is_increasing": False}, (-18, 0)),
            ({"disttype": DistType.normal, "priors": PriorClass.frequentist_unrestricted, "is_increasing": True}, (-1e6, 1e6)),
            ({"disttype": DistType.normal, "priors": PriorClass.frequentist_unrestricted, "is_increasing": False}, (-1e6, 1e6)),
            ({"disttype": DistType.normal_ncv, "priors": PriorClass.frequentist_unrestricted, "is_increasing": True}, (-18, 18)),
            ({"disttype": DistType.normal_ncv, "priors": PriorClass.frequentist_unrestricted, "is_increasing": False}, (-18, 18)),
        ]:  # fmt: on
            model = continuous.Polynomial(cdataset2, settings)
            b1 = model.settings.priors.get_prior("b1")
            assert b1.min_value == priors[0], settings
            assert b1.max_value == priors[1], settings

    def test_exp3(self, cdataset2):
        for settings, priors in [
            ({"priors": PriorClass.frequentist_restricted, "is_increasing": True}, (0, 20)),
            ({"priors": PriorClass.frequentist_restricted, "is_increasing": False}, (-20, 0)),
        ]:
            model = continuous.ExponentialM3(cdataset2, settings)
            c = model.settings.priors.get_prior("c")
            assert (c.min_value, c.max_value) == priors

    def test_exp5(self, cdataset2):
        for settings, priors in [
            ({"priors": PriorClass.frequentist_restricted, "is_increasing": True}, (0, 20)),
            ({"priors": PriorClass.frequentist_restricted, "is_increasing": False}, (-20, 0)),
        ]:
            model = continuous.ExponentialM5(cdataset2, settings)
            c = model.settings.priors.get_prior("c")
            assert (c.min_value, c.max_value) == priors

    def test_power(self, cdataset2):
        for settings, priors in [
            ({"disttype": DistType.normal}, [0.1, -100, 100]),
            ({"disttype": DistType.normal_ncv}, [1, -10000, 10000]),
        ]:
            model = continuous.Power(cdataset2, settings)
            v = model.settings.priors.get_prior("v")
            assert v.min_value == priors[1]
            assert v.max_value == priors[2]


class TestBmdModelContinuous:
    def test_get_param_names(self, cdataset2):
        # test normal model case
        for m in [
            continuous.Power(dataset=cdataset2, settings=dict(disttype=DistType.normal)),
        ]:
            assert m.get_param_names() == ["g", "v", "n", "alpha"]
        m = continuous.Power(dataset=cdataset2, settings=dict(disttype=DistType.normal_ncv))
        assert m.get_param_names() == ["g", "v", "n", "rho", "alpha"]

        # test polynomial
        model = continuous.Linear(dataset=cdataset2, settings=dict(disttype=DistType.normal))
        assert model.get_param_names() == ["g", "b1", "alpha"]
        model = continuous.Polynomial(dataset=cdataset2, settings=dict(disttype=DistType.normal))
        assert model.get_param_names() == ["g", "b1", "b2", "alpha"]
        model = continuous.Polynomial(
            dataset=cdataset2, settings=dict(degree=3, disttype=DistType.normal)
        )
        assert model.get_param_names() == ["g", "b1", "b2", "b3", "alpha"]
        model = continuous.Polynomial(
            dataset=cdataset2, settings=dict(degree=3, disttype=DistType.normal_ncv)
        )
        assert model.get_param_names() == ["g", "b1", "b2", "b3", "rho", "alpha"]

    def test_report(self, cdataset2):
        model = continuous.Hill(dataset=cdataset2)
        text = model.text()
        assert "Hill" in text
        assert "Model has not successfully executed; no results available." in text

        model.execute()
        text = model.text()
        assert "Hill" in text
        assert "Goodness of Fit:" in text

    def test_default_prior_class(self, cdataset2):
        for Model, prior_class in [
            (continuous.ExponentialM3, PriorClass.frequentist_restricted),
            (continuous.ExponentialM5, PriorClass.frequentist_restricted),
            (continuous.Power, PriorClass.frequentist_restricted),
            (continuous.Hill, PriorClass.frequentist_restricted),
            (continuous.Linear, PriorClass.frequentist_unrestricted),
            (continuous.Polynomial, PriorClass.frequentist_restricted),
        ]:
            assert Model(cdataset2).settings.priors.prior_class is prior_class

    @pytest.mark.mpl_image_compare
    def test_continuous_plot(self, cdataset2):
        model = continuous.Hill(dataset=cdataset2)
        model.execute()
        return model.plot()

    def test_automated_disttype_selection(self, cdataset, cidataset):
        # check that automated selection is correct for < 0.05
        model = continuous.Power(cdataset, settings=None)
        assert model.settings.disttype == DistType.normal_ncv
        assert model.dataset.anova().model_dump()["test2"]["TEST"] < 0.05

        # but we can override if we want
        model = continuous.Power(cdataset, settings=dict(disttype=DistType.normal))
        assert model.settings.disttype == DistType.normal

        # check that automated selection is correct for > 0.05
        model = continuous.Power(cidataset, settings=None)
        assert model.dataset.anova().model_dump()["test2"]["TEST"] > 0.05
        assert model.settings.disttype == DistType.normal

        # but we can override if we want
        model = continuous.Power(cdataset, settings=dict(disttype=DistType.normal_ncv))
        assert model.settings.disttype == DistType.normal_ncv


def test_fit_parameters(cdataset2):
    # check fit parameters for continuous modeling
    model = continuous.Power(cdataset2)
    res = model.execute()
    assert res.fit.loglikelihood == pytest.approx(res.deviance.loglikelihoods[3])


def test_increasing(cdataset2):
    """
    Basic tests to ensure AIC and BMD values are successfully created and stable for all model classes
    """
    # test increasing means dataset
    for Model, bmd_values, aic in [
        (continuous.ExponentialM3, [53, 50, 55], 3186),
        (continuous.ExponentialM5, [26, 25, 28], 3072),
        (continuous.Power, [26, 24, 30], 3070),
        (continuous.Hill, [28, 25, 32], 3072),
        (continuous.Linear, [26, 24, 28], 3068),
        (continuous.Polynomial, [26, 24, 29], 3070),
    ]:
        result = Model(cdataset2, settings=dict(disttype=DistType.normal)).execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        # for regenerating values
        # res = f"(continuous.{Model.__name__}, {np.round(actual, 0).astype(int).tolist()}, {round(result.fit.aic)}),"
        # print(res)
        assert pytest.approx(bmd_values, rel=0.05) == actual, Model.__name__
        assert pytest.approx(aic, rel=0.01) == result.fit.aic, Model.__name__


def test_decreasing(negative_cdataset):
    # test decreasing means dataset
    for Model, bmd_values, aic in [
        (continuous.ExponentialM3, [BMDS_BLANK_VALUE, BMDS_BLANK_VALUE, BMDS_BLANK_VALUE], 4298),
        (continuous.ExponentialM5, [BMDS_BLANK_VALUE, BMDS_BLANK_VALUE, BMDS_BLANK_VALUE], 4300),
        (continuous.Power, [56, 50, 64], 3080),
        (continuous.Hill, [63, 56, 70], 3086),
        (continuous.Linear, [35, 33, 38], 3117),
        (continuous.Polynomial, [52, 47, 60], 3077),
    ]:
        model = Model(negative_cdataset, settings=dict(disttype=DistType.normal))
        result = model.execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        # for regenerating values
        # res = f"(continuous.{Model.__name__}, {np.round(actual, 0).astype(int).tolist()}, {round(result.fit.aic)}),"
        # print(res)
        assert pytest.approx(bmd_values, rel=0.1) == actual, Model.__name__
        assert pytest.approx(aic, rel=0.1) == result.fit.aic, Model.__name__


def test_continuous_float_counts(cdataset2):
    # ensure float based data works
    cdataset2.ns = [n + 0.1 for n in cdataset2.ns]
    for Model, bmd_values, aic in [
        (continuous.Power, [25.9, 24.3, 29.8], 3067),
    ]:
        model = Model(cdataset2, settings=dict(disttype=DistType.normal))
        result = model.execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        assert pytest.approx(bmd_values, rel=0.05) == actual
        assert pytest.approx(aic, rel=0.01) == result.fit.aic


def test_variance(cdataset2):
    model = continuous.Power(cdataset2, dict(disttype=DistType.normal))
    result = model.execute()
    actual = [result.bmd, result.bmdl, result.bmdu]
    # print(f"{actual[0]:.1f}, {actual[1]:.1f}, {actual[2]:.1f}")
    assert model.settings.disttype is DistType.normal
    assert pytest.approx(actual, rel=0.05) == [25.9, 24.4, 29.8]
    assert result.parameters.names == ["g", "v", "n", "alpha"]

    model = continuous.Power(cdataset2, dict(disttype=DistType.normal_ncv))
    result = model.execute()
    actual = [result.bmd, result.bmdl, result.bmdu]
    # print(f"{actual[0]:.1f}, {actual[1]:.1f}, {actual[2]:.1f}")
    assert model.settings.disttype is DistType.normal_ncv
    assert result.parameters.names == ["g", "v", "n", "rho", "alpha"]
    assert pytest.approx(actual, rel=0.05) == [14.6, 13.1, 17.3]

    # only Exp can be used
    model = continuous.ExponentialM5(cdataset2, dict(disttype=DistType.log_normal))
    result = model.execute()
    actual = [result.bmd, result.bmdl, result.bmdu]
    # print(f"{actual[0]:.1f}, {actual[1]:.1f}, {actual[2]:.1f}")
    assert model.settings.disttype is DistType.log_normal
    assert pytest.approx(actual, rel=0.05) == [10.3, 9.71, 11.1]
    assert result.parameters.names == ["a", "b", "c", "d", "log-alpha"]


def test_continuous_polynomial(cdataset2):
    # compare bmd, bmdl, bmdu, aic values
    for degree, bmd_values, aic in [
        (1, [25.84, 24.358, 27.529], 3067.8),
        (2, [25.984, 24.328, 28.642], 3069.8),
        (3, [26.803, 24.259, 28.908], 3070.3),
        (4, [26.137, 24.336, 28.713], 3069.9),
    ]:
        settings = ContinuousModelSettings(degree=degree)
        result = continuous.Polynomial(cdataset2, settings).execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        # for regenerating values
        # res = f"({degree}, {np.round(actual, 3).tolist()}, {round(result.fit.aic, 1)}),"
        # print(res)
        assert pytest.approx(actual, rel=0.05) == bmd_values, degree
        assert pytest.approx(aic, rel=0.01) == result.fit.aic, degree


def test_continuous_session(cdataset2):
    session = pybmds.Session(dataset=cdataset2)
    session.add_default_models()
    session.execute()
    for model in session.models:
        model.results = model.execute()
    d = session.to_dict()
    # ensure json-serializable
    json.dumps(d)


def test_increasing_lognormal(cdataset2):
    session = pybmds.Session(dataset=cdataset2)
    settings = dict(disttype=DistType.log_normal)
    for model in (Models.ExponentialM3, Models.ExponentialM5):
        session.add_model(model, settings)
    session.execute()
    for model in session.models:
        assert model.results.has_completed is True
        assert model.results.bmd != BMDS_BLANK_VALUE

    session = pybmds.Session(dataset=cdataset2)
    settings = dict(disttype=DistType.log_normal)
    for model in (Models.Hill, Models.Power, Models.Polynomial):
        with pytest.raises(ConfigurationException):
            session.add_model(model, settings)


def test_decreasing_lognormal():
    ds = pybmds.ContinuousDataset(
        doses=[0, 10, 50, 100, 250],
        ns=[30, 30, 30, 30, 30],
        means=[0.116, 0.113, 0.108, 0.108, 0.106],
        stdevs=[0.006, 0.006, 0.004, 0.009, 0.008],
    )

    session = pybmds.Session(dataset=ds)
    settings = dict(disttype=DistType.log_normal)
    for model in (Models.ExponentialM3, Models.ExponentialM5):
        session.add_model(model, settings)
    session.execute()
    for model, bmd in zip(session.models, [227, 46], strict=True):
        assert model.results.has_completed is True
        assert model.results.bmd == pytest.approx(bmd, rel=0.05)

    session = pybmds.Session(dataset=ds)
    settings = dict(disttype=DistType.log_normal)
    for model in (Models.Hill, Models.Power, Models.Polynomial):
        with pytest.raises(ConfigurationException):
            session.add_model(model, settings)
