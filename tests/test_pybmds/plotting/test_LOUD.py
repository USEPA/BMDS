import warnings
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import xarray as xr

import pybmds
from pybmds.constants import DistType, Models, PriorClass
from pybmds.plotting.LOUD import (
    _RHAT_SINGLE_CHAIN_FOOTNOTE,
    _arviz_hdi,
    _arviz_summary,
    _as_chain_draws,
    _bmd_diagnostics_table,
    _bmd_distributions_figure,
    _bmd_summary_table,
    _build_observed_data,
    _drop_empty_summary_rows,
    _extract_model_row,
    _figure_from_axes,
    _finite_summary_table,
    _flatten_chain_draws,
    _hide_rhat_for_single_chain,
    _infer_n_chains,
    _ma_bmd_hdi,
    _ma_bmd_posterior_figure,
    _ma_bmd_quantiles,
    _model_color_map,
    _multi_summary_table,
    _pad_draws,
    _parameter_figure_height,
    _parameter_group_records,
    _parameter_group_trace_figure,
    _plot_dist,
    _rename_summary_columns,
    _reshape_draws,
    _summary_from_draws,
    get_model_average_figures,
    model_average_to_inferencedata,
)


def _data_tree(**groups: xr.Dataset) -> xr.DataTree:
    return xr.DataTree.from_dict(groups)


def _fake_loud_idata(n_chains=2):
    model_names = ["Power (CV)", "Hill (NCV)"]
    n_draws = 4
    shape = (n_chains, n_draws, len(model_names))
    coords = {"chain": np.arange(n_chains), "draw": np.arange(n_draws), "model": model_names}
    posterior = xr.Dataset(
        {
            "BMD": (
                ("chain", "draw", "model"),
                np.arange(np.prod(shape), dtype=float).reshape(shape) + 1,
            ),
            "MA_BMD": (
                ("chain", "draw"),
                np.arange(n_chains * n_draws, dtype=float).reshape(n_chains, n_draws) + 10,
            ),
            "g": (
                ("chain", "draw", "model"),
                np.arange(np.prod(shape), dtype=float).reshape(shape) + 20,
            ),
            "alpha": (
                ("chain", "draw", "model"),
                np.arange(np.prod(shape), dtype=float).reshape(shape) + 40,
            ),
            "model_weights": (("model",), np.array([0.25, 0.75])),
            "n_param": (("model",), np.array([2, 2])),
        },
        coords=coords,
    )
    return _data_tree(posterior=posterior)


def _fake_loud_session():
    class FakeModel:
        def __init__(self, name, verbose, disttype):
            self._name = name
            self.bmd_model_class = SimpleNamespace(verbose=verbose)
            self.settings = SimpleNamespace(alpha=0.05, disttype=disttype)
            self.results = SimpleNamespace(parameters=SimpleNamespace(names=["g", "alpha"]))

        def name(self):
            return self._name

    models = [
        FakeModel("Power (CV)", "Power", DistType.normal),
        FakeModel("Hill (NCV)", "Hill", DistType.normal_ncv),
    ]
    return SimpleNamespace(model_average=SimpleNamespace(models=models))


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
        idata = _data_tree(posterior=xr.Dataset(posterior, coords=coords))
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

    def test_parameter_groups_keep_lognormal_log_alpha_with_model_parameters(self):
        class FakeModel:
            bmd_model_class = SimpleNamespace(verbose="Multiplicative Hill")

            def __init__(self, name, disttype, parameter_names):
                self._name = name
                self.settings = SimpleNamespace(disttype=disttype)
                self.results = SimpleNamespace(parameters=SimpleNamespace(names=parameter_names))

            def name(self):
                return self._name

        models = [
            FakeModel("Multiplicative Hill (CV)", DistType.normal, ["a", "alpha"]),
            FakeModel(
                "Multiplicative Hill (Lognormal)",
                DistType.log_normal,
                ["a", "log-alpha"],
            ),
        ]
        session = SimpleNamespace(model_average=SimpleNamespace(models=models))
        coords = {
            "chain": [0],
            "draw": [0, 1, 2],
            "model": [model.name() for model in models],
        }
        posterior = {
            "BMD": (("chain", "draw", "model"), np.arange(6, dtype=float).reshape(1, 3, 2)),
            "a": (("chain", "draw", "model"), np.arange(10, 16, dtype=float).reshape(1, 3, 2)),
            "alpha": (
                ("chain", "draw", "model"),
                np.arange(16, 22, dtype=float).reshape(1, 3, 2),
            ),
            "log-alpha": (
                ("chain", "draw", "model"),
                np.arange(20, 26, dtype=float).reshape(1, 3, 2),
            ),
        }
        idata = _data_tree(posterior=xr.Dataset(posterior, coords=coords))

        records = _parameter_group_records(idata, session, hdi_prob=0.9)

        assert [record["name"] for record in records] == ["Multiplicative Hill"]
        main_summary = records[0]["summary"]
        assert ((main_summary["Model"] == "CV") & (main_summary["Parameter"] == "alpha")).any()
        assert (
            (main_summary["Model"] == "Lognormal") & (main_summary["Parameter"] == "log-alpha")
        ).any()

        for record in records:
            plt.close(record["trace_figure"])

    def test_parameter_groups_can_expand_by_individual_model(self):
        class FakeModel:
            bmd_model_class = SimpleNamespace(verbose="Power")

            def __init__(self, name, disttype):
                self._name = name
                self.settings = SimpleNamespace(disttype=disttype)
                self.results = SimpleNamespace(parameters=SimpleNamespace(names=["g", "alpha"]))

            def name(self):
                return self._name

        models = [
            FakeModel("Power (CV)", DistType.normal),
            FakeModel("Power (NCV)", DistType.normal_ncv),
        ]
        session = SimpleNamespace(model_average=SimpleNamespace(models=models))
        coords = {
            "chain": [0],
            "draw": [0, 1, 2],
            "model": [model.name() for model in models],
        }
        posterior = {
            "BMD": (("chain", "draw", "model"), np.arange(6, dtype=float).reshape(1, 3, 2)),
            "g": (("chain", "draw", "model"), np.arange(10, 16, dtype=float).reshape(1, 3, 2)),
            "alpha": (
                ("chain", "draw", "model"),
                np.arange(20, 26, dtype=float).reshape(1, 3, 2),
            ),
        }
        idata = _data_tree(posterior=xr.Dataset(posterior, coords=coords))

        records = _parameter_group_records(idata, session, hdi_prob=0.9, compressed=False)

        assert [record["name"] for record in records] == ["Power CV", "Power NCV"]
        for record, expected_model in zip(records, ["CV", "NCV"], strict=True):
            assert record["summary"]["Model"].unique().tolist() == [expected_model]
            assert record["var_names"] == ["g", "alpha"]
            plt.close(record["trace_figure"])

    def test_parameter_group_figures_do_not_trigger_pyplot_open_warning(self, recwarn):
        idata = _data_tree(
            posterior=xr.Dataset(
                {
                    "BMD": (
                        ("chain", "draw", "model"),
                        np.arange(3, dtype=float).reshape(1, 3, 1),
                    ),
                    "g": (
                        ("chain", "draw", "model"),
                        np.arange(3, dtype=float).reshape(1, 3, 1),
                    ),
                },
                coords={"chain": [0], "draw": [0, 1, 2], "model": ["Power (CV)"]},
            )
        )
        color_map = {"Power (CV)": "black"}

        with plt.rc_context({"figure.max_open_warning": 1}):
            figures = [
                _parameter_group_trace_figure(
                    idata=idata,
                    group_name=f"Power {idx}",
                    model_names=["Power (CV)"],
                    param_names=["g"],
                    color_map=color_map,
                )
                for idx in range(3)
            ]

        messages = [str(warning.message) for warning in recwarn]
        assert not any(
            "More than" in message and "figures have been opened" in message for message in messages
        )
        for fig in figures:
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

    def test_loud_helper_edge_branches(self, monkeypatch):
        fig, ax = plt.subplots()
        _plot_dist(np.array([1.0, 1.0]), ax=ax, color="black")
        assert ax.patches
        plt.close(fig)

        fig, ax = plt.subplots()

        def raise_kde(*args, **kwargs):
            raise ValueError("array must not contain infs or NaNs")

        monkeypatch.setattr("pybmds.plotting.LOUD.gaussian_kde", raise_kde)
        _plot_dist(np.array([1.0, 2.0, 3.0]), ax=ax, color="black")
        assert ax.patches
        plt.close(fig)

        fig, ax = plt.subplots()

        class WarningKde:
            def __init__(self, draws):
                warnings.warn("overflow encountered in dot", RuntimeWarning, stacklevel=2)

            def __call__(self, x):
                warnings.warn("overflow encountered in square", RuntimeWarning, stacklevel=2)
                return np.ones_like(x, dtype=float)

        monkeypatch.setattr("pybmds.plotting.LOUD.gaussian_kde", WarningKde)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            _plot_dist(np.array([1.0, 2.0, 3.0]), ax=ax, color="black")
        assert ax.lines
        assert not [warning for warning in caught if issubclass(warning.category, RuntimeWarning)]
        plt.close(fig)

        already_shaped = np.arange(6, dtype=float).reshape(2, 3)
        np.testing.assert_array_equal(_reshape_draws(already_shaped, n_chains=2), already_shaped)
        assert _pad_draws(np.ones((2, 1, 1)), 3).shape == (3, 1, 1)
        with pytest.raises(ValueError, match="Unsupported draw shape"):
            _pad_draws(np.ones((1, 1, 1, 1)), 2)

        np.testing.assert_array_equal(_as_chain_draws(already_shaped, 2), already_shaped)
        assert _as_chain_draws(np.arange(12, dtype=float).reshape(6, 2), 3).shape == (3, 2, 2)
        np.testing.assert_array_equal(_flatten_chain_draws(np.array([1.0, 2.0])), [1.0, 2.0])

        assert (
            _infer_n_chains(
                SimpleNamespace(model_average=SimpleNamespace(structs=SimpleNamespace(n_chains=4)))
            )
            == 4
        )
        assert (
            _infer_n_chains(
                SimpleNamespace(
                    model_average=SimpleNamespace(
                        structs=None,
                        models=[SimpleNamespace(settings=SimpleNamespace(n_chains=3))],
                    )
                )
            )
            == 3
        )
        assert (
            _infer_n_chains(SimpleNamespace(model_average=SimpleNamespace(structs=None, models=[])))
            == 1
        )

        with pytest.raises(ValueError, match="Unsupported dataset type"):
            _build_observed_data(SimpleNamespace(doses=[0, 1]))

        assert _parameter_figure_height(3) == 8.0
        assert _parameter_figure_height(6) == 11.0
        assert _parameter_figure_height(20) == 11.0

        monkeypatch.setattr("pybmds.plotting.LOUD._ARVIZ_USES_DATATREE", False)

        def fake_hdi(draws, hdi_prob=None, prob=None):
            assert hdi_prob == 0.9
            assert prob is None
            return np.array([0.0, 1.0])

        monkeypatch.setattr("pybmds.plotting.LOUD.az.hdi", fake_hdi)
        np.testing.assert_array_equal(_arviz_hdi(np.array([0.0, 1.0]), 0.9), [0.0, 1.0])

        class FakeInferenceData:
            def __init__(self, **kwargs):
                self.groups = kwargs

            @property
            def posterior(self):
                return self.groups["posterior"]

        monkeypatch.setattr("pybmds.plotting.LOUD.az.InferenceData", FakeInferenceData)
        seen = {}

        def fake_summary(data, var_names, hdi_prob=None, stat_focus=None, **kwargs):
            seen["legacy"] = (var_names, hdi_prob, stat_focus, type(data).__name__)
            return pd.DataFrame({"median": [1.0]}, index=var_names)

        monkeypatch.setattr("pybmds.plotting.LOUD.az.summary", fake_summary)
        legacy_idata = _data_tree(
            posterior=xr.Dataset(
                {"x": (("chain", "draw"), np.array([[1.0, 2.0]]))},
                coords={"chain": [0], "draw": [0, 1]},
            )
        )
        _arviz_summary(legacy_idata, ["x"], 0.9, "median")
        assert seen["legacy"] == (["x"], 0.9, "median", "FakeInferenceData")

    def test_arviz_summary_suppresses_runtime_warnings(self, monkeypatch):
        def fake_summary(*args, **kwargs):
            warnings.warn("overflow encountered in square", RuntimeWarning, stacklevel=2)
            return pd.DataFrame({"median": [1.0]}, index=["x"])

        monkeypatch.setattr("pybmds.plotting.LOUD.az.summary", fake_summary)
        idata = _data_tree(
            posterior=xr.Dataset(
                {"x": (("chain", "draw"), np.array([[1.0, 2.0]]))},
                coords={"chain": [0], "draw": [0, 1]},
            )
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            actual = _arviz_summary(idata, ["x"], 0.9, "median")

        assert actual.loc["x", "median"] == 1.0
        assert not [warning for warning in caught if issubclass(warning.category, RuntimeWarning)]

    def test_summary_helpers_empty_and_error_branches(self, monkeypatch):
        assert _summary_from_draws(np.array([np.nan]), "x", "x", 0.9).empty
        with pytest.raises(ValueError, match="Unsupported draw shape"):
            _summary_from_draws(np.ones((1, 1, 1)), "x", "x", 0.9)

        def raise_summary(*args, **kwargs):
            raise ValueError("boom")

        monkeypatch.setattr("pybmds.plotting.LOUD._arviz_summary", raise_summary)
        fallback = _summary_from_draws(np.array([[1.0, 2.0, 3.0]]), "x", "label", 0.9)
        assert fallback.loc["label", "median"] == 2.0

        idata = _data_tree(
            posterior=xr.Dataset(
                {"x": (("chain", "draw", "model", "extra"), np.ones((1, 2, 1, 1)))},
                coords={"chain": [0], "draw": [0, 1], "model": ["m"], "extra": [0]},
            )
        )
        with pytest.raises(ValueError, match="Unsupported posterior shape"):
            _finite_summary_table(idata, ["x"], 0.9)

    def test_loud_diagnostics_and_figures_with_synthetic_idata(self):
        idata = _fake_loud_idata(n_chains=2)
        session = _fake_loud_session()

        bmd_summary = _bmd_diagnostics_table(idata, hdi_prob=0.9)
        assert list(bmd_summary.index) == ["Power (CV)", "Hill (NCV)", "MA_BMD"]
        assert "BMD" in bmd_summary.columns

        parameter_groups = _parameter_group_records(idata, session, hdi_prob=0.9)
        assert [group["name"] for group in parameter_groups] == ["Power", "Hill"]
        assert parameter_groups[0]["summary"]["Model"].tolist() == ["CV", "CV"]

        stats = _ma_bmd_quantiles(idata, alpha=0.05)
        posterior = _ma_bmd_posterior_figure(idata, stats)
        overlay = _bmd_distributions_figure(idata)
        assert posterior.axes[0].get_xlabel() == "BMD"
        assert overlay.axes[0].get_title() == "BMD distributions"

        for group in parameter_groups:
            plt.close(group["trace_figure"])
        plt.close(posterior)
        plt.close(overlay)

    def test_bmd_diagnostics_falls_back_when_summary_has_no_rows(self, monkeypatch):
        idata = _fake_loud_idata(n_chains=1)
        monkeypatch.setattr(
            "pybmds.plotting.LOUD._multi_summary_table",
            lambda *args, **kwargs: pd.DataFrame(),
        )

        bmd_summary = _bmd_diagnostics_table(idata, hdi_prob=0.9)

        assert list(bmd_summary.index) == ["Power (CV)", "Hill (NCV)", "MA_BMD"]
        assert bmd_summary.loc["MA_BMD", "BMD"] == pytest.approx(11.5)

    def test_loud_figures_skip_values_too_large_for_matplotlib_layout(self):
        huge = 1e200
        model_names = ["Power"]
        idata = _data_tree(
            posterior=xr.Dataset(
                {
                    "BMD": (
                        ("chain", "draw", "model"),
                        np.array([[[huge], [huge * 1.1], [huge * 1.2]]]),
                    ),
                    "MA_BMD": (
                        ("chain", "draw"),
                        np.array([[huge, huge * 1.1, huge * 1.2]]),
                    ),
                    "g": (
                        ("chain", "draw", "model"),
                        np.array([[[huge], [huge * 1.1], [huge * 1.2]]]),
                    ),
                },
                coords={"chain": [0], "draw": [0, 1, 2], "model": model_names},
            )
        )

        posterior = _ma_bmd_posterior_figure(
            idata,
            {"lower": huge, "median": huge * 1.1, "upper": huge * 1.2},
        )
        overlay = _bmd_distributions_figure(idata)
        trace = _parameter_group_trace_figure(
            idata=idata,
            group_name="Power",
            model_names=model_names,
            param_names=["g"],
            color_map={"Power": "black"},
        )

        assert posterior.axes
        assert overlay.axes
        assert trace.axes
        plt.close(posterior)
        plt.close(overlay)
        plt.close(trace)

    def test_get_model_average_figures_hides_single_chain_rhat(self, monkeypatch):
        idata = _fake_loud_idata(n_chains=1)
        session = _fake_loud_session()

        monkeypatch.setattr("pybmds.plotting.LOUD.model_average_to_inferencedata", lambda _: idata)

        figures = get_model_average_figures(session, compressed=False)

        assert figures["bmd_summary"].attrs["footnotes"] == [_RHAT_SINGLE_CHAIN_FOOTNOTE]
        assert figures["multi_summary"].attrs["footnotes"] == [_RHAT_SINGLE_CHAIN_FOOTNOTE]
        assert "R-hat" not in figures["bmd_summary"].columns
        assert all(
            group["summary"].attrs["footnotes"] == [_RHAT_SINGLE_CHAIN_FOOTNOTE]
            for group in figures["parameter_groups"]
        )
        assert [group["name"] for group in figures["parameter_groups"]] == [
            "Power CV",
            "Hill NCV",
        ]

        for key in ("posterior", "overlay"):
            plt.close(figures[key])
        for group in figures["parameter_groups"]:
            plt.close(group["trace_figure"])

    def test_loud_row_extraction_and_rhat_footnote_helpers(self):
        summary = pd.DataFrame(
            {"r_hat": [1.0, 1.1, 1.2], "value": [1.0, 2.0, 3.0]},
            index=["x", "x[Hill]", "x[0:Power]"],
        )

        assert _extract_model_row(summary, "x", "MA_BMD")["value"] == 1.0
        assert _extract_model_row(summary, "x", "Hill")["value"] == 2.0
        assert _extract_model_row(summary, "x", "Power")["value"] == 3.0
        assert _extract_model_row(summary, "x", "Gamma") is None

        hidden = _hide_rhat_for_single_chain(summary)
        hidden = _hide_rhat_for_single_chain(hidden)
        assert "r_hat" not in hidden.columns
        assert hidden.attrs["footnotes"] == [_RHAT_SINGLE_CHAIN_FOOTNOTE]

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

    def test_build_observed_data_continuous_individual(self, cidataset):
        observed = _build_observed_data(cidataset)

        assert observed.attrs["dtype"] == "continuous_individual"
        assert observed.sizes["dose_group"] == len(cidataset.doses)
        assert observed.sizes["observation"] == len(cidataset.individual_doses)
        np.testing.assert_array_equal(observed["dose"].values, np.asarray(cidataset.doses))
        np.testing.assert_array_equal(observed["n"].values, np.asarray(cidataset.ns))
        np.testing.assert_array_equal(observed["mean"].values, np.asarray(cidataset.means))
        np.testing.assert_array_equal(observed["stdev"].values, np.asarray(cidataset.stdevs))
        np.testing.assert_array_equal(
            observed["individual_dose"].values,
            np.asarray(cidataset.individual_doses),
        )
        np.testing.assert_array_equal(observed["response"].values, np.asarray(cidataset.responses))

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

        idata = model_average_to_inferencedata(session)

        assert isinstance(idata, xr.DataTree)
        assert set(["BMD", "MA_BMD", "model_weights", "n_param"]).issubset(
            idata.posterior.data_vars
        )
        assert list(idata.posterior.coords["model"].values) == [
            model.name() for model in session.models
        ]
        assert idata.observed_data.attrs["dtype"] == "continuous"
        assert idata.posterior["MA_BMD"].shape[0] == 1
        assert idata.posterior["BMD"].shape[2] == len(session.models)

    def test_model_average_to_inferencedata_pads_filtered_json_draws(self, cdataset3):
        class FakeModel:
            bmd_model_class = SimpleNamespace(verbose="Hill")

            def __init__(self):
                self.results = SimpleNamespace(parameters=SimpleNamespace(names=["g"]))
                self.settings = SimpleNamespace(disttype=DistType.normal)

            def name(self):
                return "Hill"

        model = FakeModel()
        session = SimpleNamespace(
            dataset=cdataset3,
            model_average=SimpleNamespace(
                models=[model],
                results=SimpleNamespace(
                    bmd_dist=np.array([0.1, 0.2]),
                    model_bmd_dist=[np.array([1.0, 2.0, 3.0])],
                    model_parm_dist=[np.array([[10.0], [20.0], [30.0]])],
                    posteriors=np.array([1.0]),
                ),
            ),
        )

        idata = model_average_to_inferencedata(session)

        assert idata.posterior["MA_BMD"].shape == (1, 3)
        assert idata.posterior["BMD"].shape == (1, 3, 1)
        assert idata.posterior["g"].shape == (1, 3, 1)
        assert np.isnan(idata.posterior["MA_BMD"].values[0, 2])
        np.testing.assert_array_equal(idata.posterior["BMD"].values[0, :, 0], [1.0, 2.0, 3.0])

    def test_model_average_to_inferencedata_rejects_unpaired_model_draws(self, cdataset3):
        class FakeModel:
            bmd_model_class = SimpleNamespace(verbose="Hill")

            def __init__(self):
                self.results = SimpleNamespace(parameters=SimpleNamespace(names=["g"]))
                self.settings = SimpleNamespace(disttype=DistType.normal)

            def name(self):
                return "Hill"

        model = FakeModel()
        session = SimpleNamespace(
            dataset=cdataset3,
            model_average=SimpleNamespace(
                models=[model],
                results=SimpleNamespace(
                    bmd_dist=np.array([0.1, 0.2, 0.3]),
                    model_bmd_dist=[np.array([1.0, 2.0, 3.0])],
                    model_parm_dist=[np.array([[10.0], [20.0]])],
                    posteriors=np.array([1.0]),
                ),
            ),
        )

        with pytest.raises(ValueError, match="BMD draws but 2 parameter draws"):
            model_average_to_inferencedata(session)

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

        idata = model_average_to_inferencedata(session)

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
        idata = _data_tree(posterior=xr.Dataset(posterior, coords=coords))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            actual = _multi_summary_table(idata=idata, var_names=["a", "BMD"], hdi_prob=0.9)

        assert not [warning for warning in caught if issubclass(warning.category, RuntimeWarning)]
        assert "a[0:A]" in actual.index
        assert "a[1:B]" not in actual.index
        assert "BMD[1:B]" in actual.index

    def test_bmd_diagnostics_retains_models_with_sparse_nonfinite_draws(self):
        posterior = {
            "BMD": (
                ("chain", "draw", "model"),
                np.array(
                    [
                        [[1.0, 4.0], [2.0, np.nan], [3.0, 6.0]],
                        [[1.1, 4.1], [2.1, 5.1], [3.1, 6.1]],
                    ]
                ),
            ),
            "MA_BMD": (
                ("chain", "draw"),
                np.array([[2.0, np.nan, 4.0], [2.1, 3.1, 4.1]]),
            ),
        }
        coords = {"chain": [0, 1], "draw": np.arange(3), "model": ["A", "B"]}
        idata = _data_tree(posterior=xr.Dataset(posterior, coords=coords))

        actual = _bmd_diagnostics_table(idata, hdi_prob=0.9)

        assert list(actual.index) == ["A", "B", "MA_BMD"]
        assert actual.loc[["A", "B", "MA_BMD"], "BMD"].notna().all()

    def test_dichotomous_loud_model_average_to_inferencedata(self):
        dataset = pybmds.DichotomousDataset(
            doses=[0, 0.25, 0.75, 0.85, 1],
            ns=[20, 20, 20, 20, 20],
            incidences=[0, 1, 7, 15, 19],
        )
        session = pybmds.Session(dataset=dataset)
        session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)
        session.execute()

        idata = model_average_to_inferencedata(session)

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

        idata = model_average_to_inferencedata(session)

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
        idata = model_average_to_inferencedata(session)
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

    def test_rename_summary_columns_handles_arviz_datatree_eti_labels(self):
        summary = pd.DataFrame(
            {
                "median": [1.0],
                "eti90_lb": [0.5],
                "eti90_ub": [1.5],
            }
        )

        actual = _rename_summary_columns(summary)

        assert "Median" in actual.columns
        assert "5%" in actual.columns
        assert "95%" in actual.columns
        assert "eti90_lb" not in actual.columns
        assert "eti90_ub" not in actual.columns

    def test_rename_summary_columns_handles_arviz_datatree_eti_labels_for_bmd(self):
        summary = pd.DataFrame(
            {
                "median": [1.0],
                "eti90_lb": [0.5],
                "eti90_ub": [1.5],
            }
        )

        actual = _rename_summary_columns(summary, bmd_labels=True)

        assert "BMD" in actual.columns
        assert "BMDL" in actual.columns
        assert "BMDU" in actual.columns
        assert "eti90_lb" not in actual.columns
        assert "eti90_ub" not in actual.columns

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
        idata = model_average_to_inferencedata(session)
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

        idata = model_average_to_inferencedata(session)
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

        monkeypatch.setattr("pybmds.plotting.LOUD._arviz_summary", fake_summary)

        actual = _multi_summary_table(idata=None, var_names=["MA_BMD"], hdi_prob=0.94)

        assert "median" in actual.columns
        assert "r_hat" in actual.columns
        assert "ess_bulk" in actual.columns
        assert "ess_tail" in actual.columns
        assert "ess_median" not in actual.columns

    def test_summary_from_draws_preserves_chain_dimension(self, monkeypatch):
        seen_shapes = []

        def fake_summary(idata, var_names, hdi_prob, stat_focus):
            seen_shapes.append(idata.posterior[var_names[0]].shape)
            return pd.DataFrame(
                {
                    "median": [1.0],
                    "r_hat": [1.01],
                    "ess_bulk": [100.0],
                    "ess_tail": [90.0],
                },
                index=pd.Index([var_names[0]], name="var_name"),
            )

        monkeypatch.setattr("pybmds.plotting.LOUD._arviz_summary", fake_summary)

        actual = _summary_from_draws(
            np.arange(12, dtype=float).reshape(4, 3), "BMD", "BMD[Power]", 0.9
        )

        assert seen_shapes == [(4, 3), (4, 3)]
        assert actual.loc["BMD[Power]", "r_hat"] == 1.01

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

        monkeypatch.setattr("pybmds.plotting.LOUD._arviz_summary", fake_summary)

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
        idata = _data_tree(posterior=xr.Dataset(posterior, coords=coords))

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

        monkeypatch.setattr("pybmds.plotting.LOUD._arviz_summary", fake_summary)
        monkeypatch.setattr("pybmds.plotting.LOUD._plot_dist", fake_plot_dist)

        figures = get_model_average_figures(session)

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
        assert "R-hat" not in figures["bmd_summary"].columns
        assert figures["bmd_summary"].attrs["footnotes"] == [_RHAT_SINGLE_CHAIN_FOOTNOTE]
        assert "Bulk Effective Sample Size" in figures["bmd_summary"].columns
        assert "Tail Effective Sample Size" in figures["bmd_summary"].columns
        assert "r_hat" not in figures["multi_summary"].columns
        assert figures["multi_summary"].attrs["footnotes"] == [_RHAT_SINGLE_CHAIN_FOOTNOTE]
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
        assert {"Model", "Parameter", "Median"}.issubset(
            figures["parameter_groups"][0]["summary"].columns
        )
        assert "R-hat" not in figures["parameter_groups"][0]["summary"].columns
        assert figures["parameter_groups"][0]["summary"].attrs["footnotes"] == [
            _RHAT_SINGLE_CHAIN_FOOTNOTE
        ]
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
