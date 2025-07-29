from textwrap import dedent

import numpy as np
import pytest

import pybmds
from pybmds.constants import PriorClass
from pybmds.datasets import NestedDichotomousDataset
from pybmds.models import nested_dichotomous
from pybmds.plotting.nested_dichotomous import dose_litter_response_plot


class TestNestedLogistic:
    def test_param_names(self, nd_dataset, nd_dataset4):
        assert nd_dataset.num_dose_groups == 3
        assert nd_dataset4.num_dose_groups == 4

        model = nested_dichotomous.NestedLogistic(nd_dataset)
        names = model.get_param_names()
        assert names == ["a", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3"]

        model = nested_dichotomous.NestedLogistic(nd_dataset4)
        names = model.get_param_names()
        assert names == ["a", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3", "phi4"]

    def test_model_priors(self, nd_dataset):
        m = nested_dichotomous.NestedLogistic(dataset=nd_dataset)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤═══════╤════════╕
        │ Parameter   │   Initial │   Min │    Max │
        ╞═════════════╪═══════════╪═══════╪════════╡
        │ a           │         0 │     0 │  1     │
        │ b           │         0 │   -18 │ 18     │
        │ theta1      │         0 │     0 │  1     │
        │ theta2      │         0 │   -18 │ 18     │
        │ rho         │         1 │     1 │ 18     │
        │ phi1        │         0 │     0 │  1e+06 │
        │ phi2        │         0 │     0 │  1e+06 │
        │ phi3        │         0 │     0 │  1e+06 │
        ╘═════════════╧═══════════╧═══════╧════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check overrides
        m.settings.priors.update("a", initial_value=1, min_value=1, max_value=2)
        m.settings.priors.update("phi2", initial_value=5, min_value=5, max_value=6)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤═══════╤════════╕
        │ Parameter   │   Initial │   Min │    Max │
        ╞═════════════╪═══════════╪═══════╪════════╡
        │ a           │         1 │     1 │  2     │
        │ b           │         0 │   -18 │ 18     │
        │ theta1      │         0 │     0 │  1     │
        │ theta2      │         0 │   -18 │ 18     │
        │ rho         │         1 │     1 │ 18     │
        │ phi1        │         0 │     0 │  1e+06 │
        │ phi2        │         5 │     5 │  6     │
        │ phi3        │         0 │     0 │  1e+06 │
        ╘═════════════╧═══════════╧═══════╧════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check prior classes
        m = nested_dichotomous.NestedLogistic(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_restricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 1
        m = nested_dichotomous.NestedLogistic(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_unrestricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 0

    def test_execute(self, nd_dataset4):
        # add seed for reproducibility
        analysis = nested_dichotomous.NestedLogistic(nd_dataset4, settings=dict(bootstrap_seed=1))
        result = analysis.execute()
        assert result.has_completed is True
        assert result.bmd > 0
        text = analysis.text()
        assert len(text) > 0

    def test_float_inputs(self, nd_dataset4):
        # add floating point inputs instead of integer values
        nd_dataset4.litter_ns[0] += 0.25
        nd_dataset4.incidences[0] += 0.25
        nd_dataset4.litter_covariates[0] += 0.25
        analysis = nested_dichotomous.NestedLogistic(nd_dataset4, settings=dict(bootstrap_seed=1))
        result = analysis.execute()
        assert result.has_completed is True
        assert result.bmd > 0

    def test_penalize_aic_on_boundary(self, nd_dataset4):
        # simulate a very linear dataset
        ds = NestedDichotomousDataset(
            doses=np.repeat([1, 2, 3, 4], 5).tolist(),
            litter_ns=(np.ones(20) * 5).tolist(),
            incidences=np.repeat([1, 2, 3, 4], 5).tolist(),
            litter_covariates=np.repeat([1, 2, 3, 4], 5).tolist(),
        )
        penalized = nested_dichotomous.NestedLogistic(
            ds, settings=dict(penalize_aic_on_boundary=True)
        )
        penalized.execute()
        unpenalized = nested_dichotomous.NestedLogistic(
            ds, settings=dict(penalize_aic_on_boundary=False)
        )
        unpenalized.execute()

        assert np.isclose(
            penalized.results.parameters,
            unpenalized.results.parameters,
        ).all()
        assert np.isclose(
            penalized.results.aic - 8,  # 4 dose groups; 4 parameters penalized
            unpenalized.results.aic,
        )


class TestNctr:
    def test_param_names(self, nd_dataset, nd_dataset4):
        assert nd_dataset.num_dose_groups == 3
        assert nd_dataset4.num_dose_groups == 4

        model = nested_dichotomous.Nctr(nd_dataset)
        names = model.get_param_names()
        assert names == ["a", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3"]

        model = nested_dichotomous.Nctr(nd_dataset4)
        names = model.get_param_names()
        assert names == ["a", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3", "phi4"]

    def test_model_priors(self, nd_dataset):
        m = nested_dichotomous.Nctr(dataset=nd_dataset)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤═══════════╤═════════╕
        │ Parameter   │   Initial │       Min │     Max │
        ╞═════════════╪═══════════╪═══════════╪═════════╡
        │ a           │  0        │  0        │ 18      │
        │ b           │  0        │  0        │ 18      │
        │ theta1      │ -0.142857 │ -0.142857 │ -0.0625 │
        │ theta2      │ -0.142857 │ -0.142857 │ -0.0625 │
        │ rho         │  1        │  1        │ 18      │
        │ phi1        │  0        │  0        │  1e+08  │
        │ phi2        │  0        │  0        │  1e+08  │
        │ phi3        │  0        │  0        │  1e+08  │
        ╘═════════════╧═══════════╧═══════════╧═════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check overrides
        m.settings.priors.update("a", initial_value=1, min_value=1, max_value=2)
        m.settings.priors.update("phi2", initial_value=5, min_value=5, max_value=6)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤═══════════╤═════════╕
        │ Parameter   │   Initial │       Min │     Max │
        ╞═════════════╪═══════════╪═══════════╪═════════╡
        │ a           │  1        │  1        │  2      │
        │ b           │  0        │  0        │ 18      │
        │ theta1      │ -0.142857 │ -0.142857 │ -0.0625 │
        │ theta2      │ -0.142857 │ -0.142857 │ -0.0625 │
        │ rho         │  1        │  1        │ 18      │
        │ phi1        │  0        │  0        │  1e+08  │
        │ phi2        │  5        │  5        │  6      │
        │ phi3        │  0        │  0        │  1e+08  │
        ╘═════════════╧═══════════╧═══════════╧═════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check prior classes
        m = nested_dichotomous.Nctr(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_restricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 1
        m = nested_dichotomous.Nctr(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_unrestricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 0

    def test_execute(self, nd_dataset4):
        # add seed for reproducibility
        analysis = nested_dichotomous.Nctr(nd_dataset4, settings=dict(bootstrap_seed=1))
        result = analysis.execute()
        assert result.has_completed is True
        assert result.bmd > 0
        text = analysis.text()
        assert len(text) > 0


class TestSession:
    @pytest.mark.mpl_image_compare
    def test_dose_litter_response_plot(self, nd_dataset4):
        session = pybmds.Session(dataset=nd_dataset4)
        for settings in [
            {"litter_specific_covariate": 1, "intralitter_correlation": 1},
            {"litter_specific_covariate": 1, "intralitter_correlation": 0},
            {"litter_specific_covariate": 0, "intralitter_correlation": 1},
            {"litter_specific_covariate": 0, "intralitter_correlation": 0},
        ]:
            session.add_default_models(settings=settings)
        session.execute()
        return dose_litter_response_plot(session)
