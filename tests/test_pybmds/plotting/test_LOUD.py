import warnings
from types import SimpleNamespace

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import xarray as xr

import pybmds
from pybmds.constants import DistType, Models, PriorClass
from pybmds.plotting.LOUD import (
    _bmd_summary_table,
    _build_observed_data,
    _drop_empty_summary_rows,
    _figure_from_axes,
    _ma_bmd_hdi,
    _ma_bmd_quantiles,
    _model_color_map,
    _multi_summary_table,
    _parameter_group_trace_figure,
    _reshape_draws,
    get_model_average_figures,
    model_average_to_inferencedata,
)


class TestLOUD:
    def test_model_color_map_matches_bmds_dichotomous_legend(self):
        colors = _model_color_map(
            [
                "Hill",
                "Gamma",
                "Logistic",
                "LogLogistic",
                "LogProbit",
                "Multistage 2",
                "Probit",
                "Quantal Linear",
                "Weibull",
            ]
        )

        assert colors == {
            "Hill": "#d62728",
            "Gamma": "#1f77b4",
            "Logistic": "#2ca02c",
            "LogLogistic": "#9467bd",
            "LogProbit": "#ff7f0e",
            "Multistage 2": "#e7ba52",
            "Probit": "#8c564b",
            "Quantal Linear": "#e377c2",
            "Weibull": "#7f7f7f",
        }

    def test_parameter_trace_uses_global_model_colors(self):
        model_names = ["Power (CV)", "Hill (CV)"]
        coords = {"chain": [0], "draw": [0, 1, 2], "model": model_names}
        posterior = {
            "BMD": (("chain", "draw", "model"), np.arange(6, dtype=float).reshape(1, 3, 2)),
            "a": (("chain", "draw", "model"), np.arange(6, dtype=float).reshape(1, 3, 2)),
        }
        idata = az.InferenceData(posterior=xr.Dataset(posterior, coords=coords))
        color_map = _model_color_map(model_names)

        fig = _parameter_group_trace_figure(
            idata=idata,
            group_name="Hill",
            model_names=["Hill (CV)"],
            param_names=["a"],
            color_map=color_map,
        )

        assert fig.axes[1].lines[0].get_color() == color_map["Hill (CV)"]
        plt.close(fig)

    def test_reshape_draws(self):
        arr = np.arange(12, dtype=float).reshape(6, 2)

        actual = _reshape_draws(arr, n_chains=3)

        assert actual.shape == (3, 2, 2)
        np.testing.assert_array_equal(actual[0], arr[:2])
        np.testing.assert_array_equal(actual[2], arr[4:])

    def test_reshape_draws_requires_even_split(self):
        with pytest.raises(ValueError, match="must be divisible"):
            _reshape_draws(np.arange(10, dtype=float), n_chains=3)

    def test_build_observed_data_dichotomous(self, ddataset2):
        observed = _build_observed_data(ddataset2)

        assert observed.attrs["dtype"] == "dichotomous"
        np.testing.assert_array_equal(observed["dose"].values, np.asarray(ddataset2.doses))
        np.testing.assert_array_equal(observed["n"].values, np.asarray(ddataset2.ns))
        np.testing.assert_array_equal(
            observed["incidence"].values, np.asarray(ddataset2.incidences)
        )

    def test_build_observed_data_continuous(self, cdataset3):
        observed = _build_observed_data(cdataset3)

        assert observed.attrs["dtype"] == "continuous"
        np.testing.assert_array_equal(observed["dose"].values, np.asarray(cdataset3.doses))
        np.testing.assert_array_equal(observed["mean"].values, np.asarray(cdataset3.means))
        np.testing.assert_array_equal(observed["stdev"].values, np.asarray(cdataset3.stdevs))

    def test_model_average_to_inferencedata(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        idata = model_average_to_inferencedata(session, n_chains=1)

        assert set(["BMD", "MA_BMD", "model_weights", "n_param"]).issubset(
            idata.posterior.data_vars
        )
        assert list(idata.posterior.coords["model"].values) == [
            model.name() for model in session.models
        ]
        assert idata.observed_data.attrs["dtype"] == "continuous"
        assert idata.posterior["MA_BMD"].shape[0] == 1
        assert idata.posterior["BMD"].shape[2] == len(session.models)
        assert "alpha" in idata.posterior.data_vars
        assert "Var0" not in idata.posterior.data_vars

    def test_model_average_to_inferencedata_replaces_infinite_draws(self):
        class FakeModel:
            results = SimpleNamespace(parameters=SimpleNamespace(names=["a"]))

            def name(self):
                return "Fake"

        session = SimpleNamespace(
            dataset=SimpleNamespace(
                doses=[0, 1],
                ns=[10, 10],
                means=[1.0, 2.0],
                stdevs=[0.1, 0.2],
            ),
            model_average=SimpleNamespace(
                models=[FakeModel()],
                results=SimpleNamespace(
                    bmd_dist=np.array([1.0, np.inf, 3.0]),
                    model_bmd_dist=[np.array([1.0, -np.inf, 3.0])],
                    model_parm_dist=[np.array([[1.0], [np.inf], [3.0]])],
                    posteriors=np.array([1.0]),
                ),
            ),
        )

        idata = model_average_to_inferencedata(session, n_chains=1)

        assert np.isnan(idata.posterior["MA_BMD"].values[0, 1])
        assert np.isnan(idata.posterior["BMD"].values[0, 1, 0])
        assert np.isnan(idata.posterior["a"].values[0, 1, 0])

    def test_multi_summary_table_skips_missing_slices_without_runtime_warnings(self):
        posterior = {
            "a": (
                ("chain", "draw", "model"),
                np.array([[[1.0, np.nan], [2.0, np.nan], [3.0, np.nan]]]),
            ),
            "BMD": (
                ("chain", "draw", "model"),
                np.array([[[1.0, np.inf], [2.0, 4.0], [3.0, 5.0]]]),
            ),
        }
        coords = {"chain": [0], "draw": np.arange(3), "model": ["A", "B"]}
        idata = az.InferenceData(posterior=xr.Dataset(posterior, coords=coords))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            actual = _multi_summary_table(idata=idata, var_names=["a", "BMD"], hdi_prob=0.9)

        assert not [warning for warning in caught if issubclass(warning.category, RuntimeWarning)]
        assert "a[0:A]" in actual.index
        assert "a[1:B]" not in actual.index
        assert "BMD[1:B]" in actual.index

    def test_dichotomous_loud_model_average_to_inferencedata(self):
        dataset = pybmds.DichotomousDataset(
            doses=[0, 0.25, 0.75, 0.85, 1],
            ns=[20, 20, 20, 20, 20],
            incidences=[0, 1, 7, 15, 19],
        )
        session = pybmds.Session(dataset=dataset)
        session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)
        session.execute()

        idata = model_average_to_inferencedata(session, n_chains=1)

        assert set(["BMD", "MA_BMD", "model_weights", "n_param"]).issubset(
            idata.posterior.data_vars
        )
        assert list(idata.posterior.coords["model"].values) == [
            model.name() for model in session.models
        ]
        assert idata.observed_data.attrs["dtype"] == "dichotomous"
        assert idata.posterior["MA_BMD"].shape[0] == 1
        assert idata.posterior["BMD"].shape[2] == len(session.models)
        assert idata.posterior["BMD"].shape[1] == len(session.model_average.results.bmd_dist)

    def test_model_average_to_inferencedata_renames_nonconstant_variance_parameters(
        self, cdataset3
    ):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        idata = model_average_to_inferencedata(session, n_chains=1)

        assert "alpha" in idata.posterior.data_vars
        assert "rho" in idata.posterior.data_vars
        assert "Var0" not in idata.posterior.data_vars
        assert "Var1" not in idata.posterior.data_vars

    def test_loud_results_use_raw_continuous_parameter_names(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.execute()

        assert session.models[0].results.parameters.names == ["g", "v", "n", "rho", "alpha"]
        assert session.models[1].results.parameters.names == ["g", "v", "k", "n", "rho", "alpha"]

    def test_bmd_summary_table_uses_bmd_labels(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        alpha = session.models[0].settings.alpha
        idata = model_average_to_inferencedata(session, n_chains=1)
        summary = _bmd_summary_table(idata, alpha)

        assert "BMD" in summary.columns
        assert "BMDL" in summary.columns
        assert "BMDU" in summary.columns
        assert "R-hat" in summary.columns
        assert "Bulk Effective Sample Size" in summary.columns
        assert "Tail Effective Sample Size" in summary.columns
        assert "median" not in summary.columns
        assert "mean" not in summary.columns
        assert "MA_BMD" in summary.index

    def test_ma_bmd_quantiles(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        alpha = session.models[0].settings.alpha
        idata = model_average_to_inferencedata(session, n_chains=1)
        quantiles = _ma_bmd_quantiles(idata, alpha)

        assert set(quantiles) == {"lower", "median", "upper"}
        assert quantiles["lower"] <= quantiles["median"] <= quantiles["upper"]

    def test_ma_bmd_hdi(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        idata = model_average_to_inferencedata(session, n_chains=1)
        stats = _ma_bmd_hdi(idata, 0.9)

        assert set(stats) == {"lower", "median", "upper"}
        assert stats["lower"] <= stats["median"] <= stats["upper"]

    def test_multi_summary_table_keeps_ess_bulk_and_ess_tail(self, monkeypatch):
        expected_index = pd.Index(["MA_BMD"], name="var_name")

        def fake_summary(idata, var_names, hdi_prob, stat_focus=None):
            if stat_focus == "median":
                return pd.DataFrame(
                    {
                        "median": [1.2],
                        "mad": [0.1],
                        "eti_3%": [1.0],
                        "eti_97%": [1.4],
                        "mcse_median": [0.01],
                        "ess_median": [100.0],
                        "r_hat": [1.0],
                    },
                    index=expected_index,
                )
            return pd.DataFrame(
                {
                    "mean": [1.2],
                    "sd": [0.1],
                    "ess_bulk": [90.0],
                    "ess_tail": [80.0],
                },
                index=expected_index,
            )

        monkeypatch.setattr("pybmds.plotting.LOUD.az.summary", fake_summary)

        actual = _multi_summary_table(idata=None, var_names=["MA_BMD"], hdi_prob=0.94)

        assert "median" in actual.columns
        assert "r_hat" in actual.columns
        assert "ess_bulk" in actual.columns
        assert "ess_tail" in actual.columns
        assert "ess_median" not in actual.columns

    def test_multi_summary_table_summarizes_each_variable_independently(self, monkeypatch):
        calls = []

        def fake_summary(idata, var_names, hdi_prob, stat_focus=None):
            calls.append((tuple(var_names), stat_focus))
            return pd.DataFrame(
                {
                    "median" if stat_focus == "median" else "mean": [1.0],
                    "ess_bulk": [90.0],
                    "ess_tail": [80.0],
                },
                index=pd.Index([var_names[0]], name="var_name"),
            )

        monkeypatch.setattr("pybmds.plotting.LOUD.az.summary", fake_summary)

        actual = _multi_summary_table(
            idata=None, var_names=["BMD", "MA_BMD", "alpha"], hdi_prob=0.9
        )

        assert calls == [
            (("BMD",), "median"),
            (("BMD",), "mean"),
            (("MA_BMD",), "median"),
            (("MA_BMD",), "mean"),
            (("alpha",), "median"),
            (("alpha",), "mean"),
        ]
        assert list(actual.index) == ["BMD", "MA_BMD", "alpha"]

    def test_multi_summary_table_falls_back_when_arviz_summary_raises(self):
        posterior = {
            "BMD": (("chain", "draw", "model"), np.arange(12, dtype=float).reshape(1, 4, 3)),
            "MA_BMD": (("chain", "draw"), np.arange(4, dtype=float).reshape(1, 4)),
        }
        coords = {"chain": [0], "draw": np.arange(4), "model": ["A", "B", "C"]}
        idata = az.InferenceData(posterior=xr.Dataset(posterior, coords=coords))

        actual = _multi_summary_table(idata=idata, var_names=["BMD", "MA_BMD"], hdi_prob=0.9)

        assert isinstance(actual, pd.DataFrame)
        assert "median" in actual.columns
        assert "r_hat" in actual.columns
        assert "ess_bulk" in actual.columns
        assert "ess_tail" in actual.columns

    def test_drop_empty_summary_rows_removes_fully_blank_rows(self):
        df = pd.DataFrame(
            {
                "median": [1.0, np.nan, 2.0],
                "mad": [0.1, np.nan, 0.2],
                "ess_bulk": [10.0, np.nan, 20.0],
            },
            index=["kept_a", "drop_me", "kept_b"],
        )

        actual = _drop_empty_summary_rows(df)

        assert list(actual.index) == ["kept_a", "kept_b"]

    def test_get_model_average_figures(self, cdataset3, monkeypatch):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        calls = {"summary": [], "dist": 0}

        def fake_summary(idata, var_names, hdi_prob, stat_focus=None):
            calls["summary"].append((tuple(var_names), hdi_prob, stat_focus))
            var_name = var_names[0]
            if var_name == "BMD":
                index = pd.Index(
                    [f"BMD[{idx}:{model.name()}]" for idx, model in enumerate(session.models)],
                    name="var_name",
                )
            elif var_name == "MA_BMD":
                index = pd.Index(["MA_BMD"], name="var_name")
            else:
                index = pd.Index(
                    [f"{var_name}[0:{session.models[0].name()}]"],
                    name="var_name",
                )
            if stat_focus == "median":
                return pd.DataFrame(
                    {
                        "median": [1.0] * len(index),
                        "mad": [0.1] * len(index),
                        "eti_5%": [0.9] * len(index),
                        "eti_95%": [1.1] * len(index),
                        "mcse_median": [0.01] * len(index),
                        "ess_median": [100.0] * len(index),
                        "r_hat": [1.0] * len(index),
                    },
                    index=index,
                )
            return pd.DataFrame(
                {
                    "mean": [1.0] * len(index),
                    "sd": [0.1] * len(index),
                    "ess_bulk": [90.0] * len(index),
                    "ess_tail": [80.0] * len(index),
                },
                index=index,
            )

        def fake_plot_dist(*args, **kwargs):
            calls["dist"] += 1
            ax = kwargs["ax"]
            color = kwargs.get("color", "black")
            ax.plot([0, 1], [0, 1], color=color)
            return ax

        monkeypatch.setattr("pybmds.plotting.LOUD.az.summary", fake_summary)
        monkeypatch.setattr("pybmds.plotting.LOUD.az.plot_dist", fake_plot_dist)

        figures = get_model_average_figures(session, n_chains=1)

        alpha = session.models[0].settings.alpha
        assert figures["alpha"] == pytest.approx(alpha)
        assert figures["hdi_prob"] == pytest.approx(0.9)
        assert set(figures["ma_bmd_quantiles"]) == {"lower", "median", "upper"}
        assert set(figures["ma_bmd_hdi"]) == {"lower", "median", "upper"}
        assert isinstance(figures["bmd_summary"], pd.DataFrame)
        assert isinstance(figures["multi_summary"], pd.DataFrame)
        assert isinstance(figures["parameter_groups"], list)
        assert "BMD" in figures["bmd_summary"].columns
        assert "BMDL" in figures["bmd_summary"].columns
        assert "BMDU" in figures["bmd_summary"].columns
        assert "median" not in figures["bmd_summary"].columns
        assert "eti_5%" not in figures["bmd_summary"].columns
        assert "eti_95%" not in figures["bmd_summary"].columns
        assert "mean" not in figures["bmd_summary"].columns
        assert "R-hat" in figures["bmd_summary"].columns
        assert "Bulk Effective Sample Size" in figures["bmd_summary"].columns
        assert "Tail Effective Sample Size" in figures["bmd_summary"].columns
        assert "ess_bulk" in figures["multi_summary"].columns
        assert "ess_tail" in figures["multi_summary"].columns
        assert "ess_median" not in figures["multi_summary"].columns
        expected_summary_calls = []
        expected_var_names = ["BMD", "MA_BMD"]
        for model in session.models:
            for name in model.results.parameters.names:
                if name not in expected_var_names:
                    expected_var_names.append(name)
        for name in expected_var_names:
            expected_summary_calls.append(((name,), 0.9, "median"))
            expected_summary_calls.append(((name,), 0.9, "mean"))
        for expected_call in expected_summary_calls:
            assert expected_call in calls["summary"]
        assert calls["dist"] > 0
        posterior_ax = figures["posterior"].axes[0]
        legend_labels = [text.get_text() for text in posterior_ax.get_legend().get_texts()]
        assert legend_labels == ["BMDL", "BMDU", "BMD"]
        assert len(posterior_ax.lines) == 4
        annotation_text = [text.get_text() for text in posterior_ax.texts]
        assert any(text.startswith("BMDL:") for text in annotation_text)
        assert any(text.startswith("BMD:") for text in annotation_text)
        assert any(text.startswith("BMDU:") for text in annotation_text)
        assert isinstance(figures["overlay"], plt.Figure)
        assert figures["parameter_groups"][0]["name"] == "Power"
        assert {"Model", "Parameter", "Median", "R-hat"}.issubset(
            figures["parameter_groups"][0]["summary"].columns
        )
        assert len(figures["parameter_groups"][0]["trace_figure"].axes[0].lines) >= 1
        assert figures["overlay"].legends
        assert figures["parameter_groups"][0]["trace_figure"].legends
        parameter_summary = figures["parameter_groups"][0]["summary"]
        assert "Median Absolute Deviation" in parameter_summary.columns
        assert "Markov Chain Standard Error (Median)" in parameter_summary.columns
        assert "Bulk Effective Sample Size" in parameter_summary.columns
        assert "Tail Effective Sample Size" in parameter_summary.columns
        assert "5%" in parameter_summary.columns
        assert "95%" in parameter_summary.columns
        assert "mad" not in parameter_summary.columns
        assert "mcse_median" not in parameter_summary.columns
        assert "eti_5%" not in parameter_summary.columns

        plt.close(figures["posterior"])
        plt.close(figures["overlay"])
        for group in figures["parameter_groups"]:
            plt.close(group["trace_figure"])

    def test_figure_from_axes(self):
        fig, ax = plt.subplots()

        actual = _figure_from_axes(np.array([[ax]]))

        assert actual is fig
        plt.close(fig)
