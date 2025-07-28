import json

import numpy as np
import pytest

import pybmds
from pybmds.constants import BMDS_BLANK_VALUE, PriorClass, PriorDistribution
from pybmds.models import dichotomous
from pybmds.models.dichotomous import Multistage, QuantalLinear
from pybmds.models.multi_tumor import MultistageCancer
from pybmds.types.dichotomous import DichotomousModelSettings, DichotomousRiskType
from pybmds.types.priors import ModelPriors, Prior


class TestBmdModelDichotomous:
    def test_get_param_names(self, ddataset2):
        # test normal model case
        model = dichotomous.Gamma(dataset=ddataset2)
        assert model.get_param_names() == ["g", "a", "b"]

        # test multistage
        model = dichotomous.Multistage(dataset=ddataset2)
        assert model.get_param_names() == ["g", "b1", "b2"]
        model = dichotomous.Multistage(dataset=ddataset2, settings=dict(degree=3))
        assert model.get_param_names() == ["g", "b1", "b2", "b3"]

    def test_default_prior_class(self, ddataset2):
        for Model, prior_class in [
            # restricted
            (dichotomous.DichotomousHill, PriorClass.frequentist_restricted),
            (dichotomous.Gamma, PriorClass.frequentist_restricted),
            (dichotomous.LogLogistic, PriorClass.frequentist_restricted),
            (dichotomous.LogLogistic, PriorClass.frequentist_restricted),
            (dichotomous.Multistage, PriorClass.frequentist_restricted),
            (dichotomous.Weibull, PriorClass.frequentist_restricted),
            # unrestricted
            (dichotomous.Logistic, PriorClass.frequentist_unrestricted),
            (dichotomous.LogProbit, PriorClass.frequentist_unrestricted),
            (dichotomous.Probit, PriorClass.frequentist_unrestricted),
            (dichotomous.QuantalLinear, PriorClass.frequentist_unrestricted),
        ]:
            assert Model(ddataset2).settings.priors.prior_class is prior_class

    def test_report(self, ddataset2):
        model = dichotomous.Gamma(dataset=ddataset2)
        text = model.text()
        assert "Gamma" in text
        assert "Model has not successfully executed; no results available." in text

        model.execute()
        text = model.text()
        assert "Gamma" in text
        assert "Goodness of Fit:" in text

    def test_risk_type(self, ddataset2):
        # extra (default)
        model = dichotomous.Logistic(dataset=ddataset2)
        resp1 = model.execute()
        assert model.settings.bmr_type is DichotomousRiskType.ExtraRisk

        # added
        model = dichotomous.Logistic(dataset=ddataset2, settings=dict(bmr_type=0))
        resp2 = model.execute()
        assert model.settings.bmr_type is DichotomousRiskType.AddedRisk

        assert not np.isclose(resp1.bmd, resp2.bmd)
        assert resp1.bmd < resp2.bmd

    @pytest.mark.mpl_image_compare
    def test_dichotomous_plot(self, ddataset2):
        model = dichotomous.Logistic(dataset=ddataset2)
        model.execute()
        return model.plot(figsize=(5, 4))

    @pytest.mark.mpl_image_compare
    def test_dichotomous_plot_axlines(self, ddataset2):
        model = dichotomous.Logistic(dataset=ddataset2)
        model.execute()
        return model.plot(figsize=(5, 4), axlines=True)

    @pytest.mark.mpl_image_compare
    def test_dichotomous_cdf_plot(self, ddataset2):
        model = dichotomous.Logistic(dataset=ddataset2)
        model.execute()
        return model.cdf_plot()

    def test_penalize_aic_on_boundary(self, ddataset):
        penalized = Multistage(ddataset, settings=dict(penalize_aic_on_boundary=True))
        penalized.execute()
        unpenalized = Multistage(ddataset, settings=dict(penalize_aic_on_boundary=False))
        unpenalized.execute()

        assert np.isclose(
            penalized.results.parameters.values,
            unpenalized.results.parameters.values,
        ).all()
        assert np.isclose(
            penalized.results.fit.aic + 2,
            unpenalized.results.fit.aic,
        )


def test_dichotomous_models(ddataset2):
    # compare bmd, bmdl, bmdu, aic values
    for Model, bmd_values, aic in [
        (dichotomous.Logistic, [69.584, 61.193, 77.945], 364.0),
        (dichotomous.LogLogistic, [68.361, 59.795, 76.012], 365.0),
        (dichotomous.Probit, [66.883, 58.333, 75.368], 362.1),
        (dichotomous.LogProbit, [66.149, 58.684, 73.16], 366.3),
        (dichotomous.Gamma, [66.061, 57.639, 73.716], 361.6),
        (dichotomous.QuantalLinear, [17.679, 15.645, 20.062], 425.6),
        (dichotomous.Weibull, [64.26, 55.219, 72.815], 358.4),
        (dichotomous.DichotomousHill, [68.178, 59.795, 75.999], 364.0),
    ]:
        model = Model(ddataset2)
        result = model.execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        # for regenerating values
        # print(
        #     f"(dichotomous.{Model.__name__}, {np.round(actual, 3).tolist()}, {round(result.fit.aic, 1)})"
        # )
        assert pytest.approx(bmd_values, abs=0.5) == actual
        assert pytest.approx(aic, abs=3.0) == result.fit.aic


def test_dichotomous_float_counts(ddataset3):
    # ensure float based data works
    for Model, bmd_values, aic in [
        (dichotomous.Logistic, [47.2, 30.7, 69.5], 49.2),
    ]:
        model = Model(ddataset3)
        result = model.execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        assert pytest.approx(bmd_values, abs=0.5) == actual
        assert pytest.approx(aic, abs=3.0) == result.fit.aic


def test_dichotomous_multistage(ddataset2):
    # compare bmd, bmdl, bmdu, aic values
    for degree, bmd_values, aic in [
        (1, [17.680, 15.645, 20.062], 425.6),
        (2, [48.016, 44.136, 51.240], 369.7),
        (3, [63.873, 52.260, 72.126], 358.5),
        (4, [63.871, 52.073, 72.725], 358.5),
    ]:
        settings = DichotomousModelSettings(degree=degree)
        model = dichotomous.Multistage(ddataset2, settings)
        result = model.execute()
        actual = [result.bmd, result.bmdl, result.bmdu]
        # for modifying values
        # print(f"({degree}, {np.round(actual, 3).tolist()}, {round(result.fit.aic, 1)})")
        assert pytest.approx(bmd_values, abs=0.5) == actual
        assert pytest.approx(aic, abs=5.0) == result.fit.aic


def test_slope_factor(ddataset2):
    model = dichotomous.Multistage(ddataset2, settings=dict(degree=2))
    model.execute()
    assert model.results.slope_factor == BMDS_BLANK_VALUE
    assert "Slope Factor" not in model.text()

    model = dichotomous.Multistage(ddataset2, settings=dict(degree=2))
    model.execute(slope_factor=True)
    assert model.results.slope_factor > 0
    assert "Slope Factor" in model.text()


def test_dichotomous_session(ddataset2):
    session = pybmds.Session(dataset=ddataset2)
    session.add_default_models()
    session.execute()
    d = session.to_dict()
    # ensure json-serializable
    json.dumps(d)


def test_dichotomous_fit_parameters(ddataset2):
    # check fit parameters for dichotomous modeling
    model = dichotomous.Logistic(ddataset2)
    res = model.execute()
    # overall fit
    actual = [res.fit.loglikelihood, res.fit.aic, res.gof.p_value, res.gof.df, res.fit.chisq]
    assert actual == pytest.approx([-179.98, 363.96, 0.48, 3.0, 2.45], abs=0.01)
    # scaled residuals
    assert res.gof.residual == pytest.approx([-1.08, -0.42, 0.94, -0.14, -0.46], abs=0.01)
    # deviance
    assert res.deviance.deviance == pytest.approx([-9999.0, 3.57, 307.681], abs=0.01)
    assert res.deviance.df == pytest.approx([-9999, 3, 4])
    assert res.deviance.p_value == pytest.approx([-9999.0, 0.311, 0.0], abs=0.01)
    assert res.fit.loglikelihood == pytest.approx(res.deviance.ll[1])


def test_dichotomous_pvalue():
    ds = pybmds.datasets.DichotomousDataset(
        doses=[0, 10, 30, 100], ns=[20, 20, 20, 20], incidences=[0, 0, 8, 20]
    )
    m = dichotomous.Logistic(dataset=ds)
    m.execute()

    # fix case found in 2023.03 where if p_value is exactly one, would incorrectly return -9999
    assert m.results.gof.p_value == pytest.approx(1.0, abs=1e-3)


class TestMultistageCancer:
    def test_settings(self, ddataset2):
        default = json.loads(
            '{"prior_class": 1, "priors": [{"name": "g", "type": 0, "initial_value": -17.0, "stdev": 0.0, "min_value": -18.0, "max_value": 18.0}, {"name": "b1", "type": 0, "initial_value": 0.1, "stdev": 0.0, "min_value": 0.0, "max_value": 10000.0}, {"name": "bN", "type": 0, "initial_value": 0.1, "stdev": 0.0, "min_value": 0.0, "max_value": 10000.0}], "variance_priors": null, "overrides": null}'
        )

        # default MultistageCancer use cancer prior
        model = MultistageCancer(ddataset2)
        assert model.settings.bmr == 0.1
        assert model.settings.priors.model_dump() == default

        # MultistageCancer w/ custom settings, but unspecified prior use cancer prior
        model = MultistageCancer(ddataset2, settings=dict(bmr=0.2))
        assert model.settings.bmr == 0.2
        assert model.settings.priors.model_dump() == default

        # MultistageCancer w/ DichotomousModelSettings is preserved
        base_model = Multistage(ddataset2)
        model = MultistageCancer(ddataset2, DichotomousModelSettings())
        assert model.settings.bmr == 0.1
        assert model.settings.priors.model_dump() != default
        assert model.settings.priors.model_dump() == base_model.settings.priors.model_dump()

        # MultistageCancer w/ custom settings is preserved
        custom = ModelPriors(
            prior_class=PriorClass.frequentist_restricted,
            priors=[
                Prior(
                    name="g",
                    type=PriorDistribution.Uniform,
                    initial_value=0,
                    stdev=0,
                    min_value=-1,
                    max_value=1,
                )
            ],
            variance_priors=None,
        )
        model = MultistageCancer(ddataset2, settings=dict(bmr=0.2, priors=custom))
        assert model.settings.bmr == 0.2
        assert model.settings.priors.model_dump() == custom.model_dump()


class TestQuantalLinear:
    def test_dr_curve(self, ddataset):
        # check that the background param if nonzero has a y value > 0 at zero dose
        model = QuantalLinear(dataset=ddataset)
        ys = model.dr_curve(doses=np.array([0]), params=(0.15, 0.2))
        assert np.allclose(ys, [0.15])
