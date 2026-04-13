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
    _multi_summary_table,
    _reshape_draws,
    get_model_average_figures,
    model_average_to_inferencedata,
)


class TestLOUD:
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

        assert session.models[0].results.parameters.names == ["g", "v", "n", "alpha", "rho"]
        assert session.models[1].results.parameters.names == ["g", "v", "k", "n", "alpha", "rho"]

    def test_bmd_summary_table_uses_median_and_alpha_bounds(self, cdataset3):
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

        assert "median" in summary.columns
        assert f"alpha_lower ({alpha:.2f})" in summary.columns
        assert f"alpha_upper ({1 - alpha:.2f})" in summary.columns
        assert "r_hat" in summary.columns
        assert "ess_bulk" in summary.columns
        assert "ess_tail" in summary.columns
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

        calls = {"summary": [], "posterior": [], "dist": 0}

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

        def fake_plot_posterior(idata, var_names):
            calls["posterior"].append(tuple(var_names))
            fig, ax = plt.subplots()
            ax.set_xlabel("BMD")
            ax.plot([0, 1], [0, 1], label="density")
            return np.array([ax])

        def fake_plot_dist(*args, **kwargs):
            calls["dist"] += 1
            ax = kwargs["ax"]
            color = kwargs.get("color", "black")
            ax.plot([0, 1], [0, 1], color=color)
            return ax

        monkeypatch.setattr("pybmds.plotting.LOUD.az.summary", fake_summary)
        monkeypatch.setattr("pybmds.plotting.LOUD.az.plot_posterior", fake_plot_posterior)
        monkeypatch.setattr("pybmds.plotting.LOUD.az.plot_dist", fake_plot_dist)

        figures = get_model_average_figures(session, n_chains=1)

        alpha = session.models[0].settings.alpha
        lower_col = "eti_5%"
        upper_col = "eti_95%"

        assert figures["alpha"] == pytest.approx(alpha)
        assert figures["hdi_prob"] == pytest.approx(0.9)
        assert set(figures["ma_bmd_quantiles"]) == {"lower", "median", "upper"}
        assert set(figures["ma_bmd_hdi"]) == {"lower", "median", "upper"}
        assert isinstance(figures["bmd_summary"], pd.DataFrame)
        assert isinstance(figures["multi_summary"], pd.DataFrame)
        assert isinstance(figures["parameter_groups"], list)
        assert "median" in figures["bmd_summary"].columns
        assert lower_col in figures["bmd_summary"].columns
        assert upper_col in figures["bmd_summary"].columns
        assert "mean" not in figures["bmd_summary"].columns
        assert "r_hat" in figures["bmd_summary"].columns
        assert "ess_bulk" in figures["multi_summary"].columns
        assert "ess_tail" in figures["multi_summary"].columns
        assert "ess_median" not in figures["multi_summary"].columns
        assert calls["posterior"] == [("MA_BMD",)]
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
        assert len(figures["posterior"].axes[0].lines) >= 4
        assert isinstance(figures["overlay"], plt.Figure)
        assert figures["parameter_groups"][0]["name"] == "Power"
        assert {"Model", "Parameter", "median", "r_hat", "ess_bulk", "ess_tail"}.issubset(
            figures["parameter_groups"][0]["summary"]
        )
        assert len(figures["parameter_groups"][0]["trace_figure"].axes[0].lines) >= 1
        assert figures["overlay"].legends
        assert figures["parameter_groups"][0]["trace_figure"].legends

        plt.close(figures["posterior"])
        plt.close(figures["overlay"])
        for group in figures["parameter_groups"]:
            plt.close(group["trace_figure"])

    def test_figure_from_axes(self):
        fig, ax = plt.subplots()

        actual = _figure_from_axes(np.array([[ax]]))

        assert actual is fig
        plt.close(fig)
