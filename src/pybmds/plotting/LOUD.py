from __future__ import annotations

from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D

from ..constants import DistType


def _reshape_draws(arr: np.ndarray, n_chains: int) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = arr.shape[0]
    if n % n_chains != 0:
        raise ValueError(f"{n=} must be divisible by {n_chains=}.")
    draws_per_chain = n // n_chains
    return arr.reshape((n_chains, draws_per_chain) + arr.shape[1:])


def _build_observed_data(dataset) -> xr.Dataset:
    """
    Build ArviZ observed_data from a pybmds dataset.

    Supports both dichotomous and continuous grouped datasets.
    Produces a structure like:

    /observed_data
      Dimensions: (dose_group: N)
      Coordinates:
        dose_group
      Data variables:
        dose
        ...
      Attributes:
        dtype: dichotomous|continuous
    """
    if not hasattr(dataset, "doses"):
        raise ValueError("session.dataset must have a 'doses' attribute.")

    doses = np.asarray(dataset.doses, dtype=float)
    n_groups = len(doses)

    observed = xr.Dataset(
        coords={
            "dose_group": np.arange(n_groups, dtype=int),
        }
    )

    observed["dose"] = ("dose_group", doses)

    # Dichotomous dataset
    if hasattr(dataset, "incidences") and hasattr(dataset, "ns"):
        observed["n"] = ("dose_group", np.asarray(dataset.ns, dtype=int))
        observed["incidence"] = ("dose_group", np.asarray(dataset.incidences, dtype=int))
        observed.attrs["dtype"] = "dichotomous"
        return observed

    # Continuous grouped dataset
    if hasattr(dataset, "ns") and hasattr(dataset, "means") and hasattr(dataset, "stdevs"):
        observed["n"] = ("dose_group", np.asarray(dataset.ns, dtype=int))
        observed["mean"] = ("dose_group", np.asarray(dataset.means, dtype=float))
        observed["stdev"] = ("dose_group", np.asarray(dataset.stdevs, dtype=float))
        observed.attrs["dtype"] = "continuous"
        return observed

    raise ValueError(
        "Unsupported dataset type. Expected either:\n"
        "- dichotomous: doses, ns, incidences\n"
        "- continuous: doses, ns, means, stdevs"
    )


def model_average_to_inferencedata(
    session,
    n_chains: int = 1,
    path: str | Path | None = None,
) -> az.InferenceData:
    """
    Convert a pybmds session with continuous model averaging results into an
    ArviZ InferenceData object.

    Posterior group:
        - one variable per pybmds parameter name as returned by model results
        - BMD           : model-specific BMD draws
        - MA_BMD        : model-averaged BMD draws
        - model_weights : posterior model weights
        - n_param       : number of parameters per model

    Observed data group:
        - built from session.dataset
    """
    ma_result = session.model_average.results
    models = session.model_average.models

    model_names = [model.name() for model in models]
    param_names_by_model = [list(model.results.parameters.names) for model in models]

    all_param_names: list[str] = []
    for names in param_names_by_model:
        for name in names:
            if name not in all_param_names:
                all_param_names.append(name)

    n_draws = len(ma_result.bmd_dist)
    n_models = len(model_names)
    n_params = len(all_param_names)

    model_bmd = np.full((n_draws, n_models), np.nan, dtype=float)
    model_params = np.full((n_draws, n_models, n_params), np.nan, dtype=float)
    n_param_by_model = np.zeros(n_models, dtype=int)

    for model_idx, (bmds, parms, pnames) in enumerate(
        zip(ma_result.model_bmd_dist, ma_result.model_parm_dist, param_names_by_model, strict=True)
    ):
        bmds = np.asarray(bmds, dtype=float)
        parms = np.asarray(parms, dtype=float)

        if len(bmds) != n_draws:
            raise ValueError(
                f"Model {model_names[model_idx]!r} has {len(bmds)} BMD draws, expected {n_draws}."
            )

        if parms.ndim != 2:
            raise ValueError(
                f"Model {model_names[model_idx]!r} parameter draws must be 2D; got shape {parms.shape}."
            )

        if parms.shape[0] != n_draws:
            raise ValueError(
                f"Model {model_names[model_idx]!r} has {parms.shape[0]} parameter draws, expected {n_draws}."
            )

        if parms.shape[1] != len(pnames):
            raise ValueError(
                f"Model {model_names[model_idx]!r} has {parms.shape[1]} parameters but "
                f"{len(pnames)} names."
            )

        model_bmd[:, model_idx] = bmds
        n_param_by_model[model_idx] = len(pnames)

        for param_idx, pname in enumerate(pnames):
            global_idx = all_param_names.index(pname)
            model_params[:, model_idx, global_idx] = parms[:, param_idx]

    posterior_raw = xr.Dataset(
        data_vars={
            "params": (
                ("chain", "draw", "model", "param"),
                _reshape_draws(model_params, n_chains),
            ),
            "BMD": (
                ("chain", "draw", "model"),
                _reshape_draws(model_bmd, n_chains),
            ),
            "MA_BMD": (
                ("chain", "draw"),
                _reshape_draws(np.asarray(ma_result.bmd_dist, dtype=float), n_chains),
            ),
            "weights": (
                ("model",),
                np.asarray(ma_result.posteriors, dtype=float),
            ),
            "n_param": (
                ("model",),
                n_param_by_model,
            ),
        },
        coords={
            "chain": np.arange(n_chains),
            "draw": np.arange(n_draws // n_chains),
            "model": model_names,
            "param": all_param_names,
        },
    )

    model_coord_posterior = xr.Dataset(
        {
            **{
                pname: posterior_raw["params"].sel(param=pname).drop_vars("param")
                for pname in all_param_names
            },
            "BMD": posterior_raw["BMD"],
            "MA_BMD": posterior_raw["MA_BMD"],
            "model_weights": posterior_raw["weights"],
            "n_param": posterior_raw["n_param"],
        },
        coords={
            "chain": posterior_raw.chain,
            "draw": posterior_raw.draw,
            "model": posterior_raw.model,
        },
    )

    observed_data = _build_observed_data(session.dataset)

    idata = az.InferenceData(
        posterior=model_coord_posterior,
        observed_data=observed_data,
    )

    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        idata.to_netcdf(path)

    return idata


def _figure_from_axes(axes):
    """Extract a matplotlib Figure from ArviZ axes output."""
    import numpy as np

    ax = np.ravel(np.atleast_1d(axes))[0]
    return ax.figure


def _ma_bmd_quantiles(idata: az.InferenceData, alpha: float) -> dict[str, float]:
    ma_bmd = np.asarray(idata.posterior["MA_BMD"].values, dtype=float).reshape(-1)
    return {
        "lower": float(np.quantile(ma_bmd, alpha)),
        "median": float(np.quantile(ma_bmd, 0.5)),
        "upper": float(np.quantile(ma_bmd, 1 - alpha)),
    }


def _ma_bmd_hdi(idata: az.InferenceData, hdi_prob: float) -> dict[str, float]:
    ma_bmd = np.asarray(idata.posterior["MA_BMD"].values, dtype=float).reshape(-1)
    lower, upper = np.asarray(az.hdi(ma_bmd, hdi_prob=hdi_prob), dtype=float)
    return {
        "lower": float(lower),
        "median": float(np.quantile(ma_bmd, 0.5)),
        "upper": float(upper),
    }


def _add_ma_bmd_reference_lines(
    fig: plt.Figure,
    stats: dict[str, float],
    lower_label: str,
    upper_label: str,
):
    labels = {
        "lower": lower_label,
        "median": "MA_BMD median",
        "upper": upper_label,
    }
    colors = {
        "lower": "#B04A3A",
        "median": "#111111",
        "upper": "#2A7F62",
    }
    linestyles = {
        "lower": "--",
        "median": "-",
        "upper": "--",
    }

    for ax in fig.axes:
        if not ax.has_data():
            continue

        xlabel = (ax.get_xlabel() or "").lower()
        ylabel = (ax.get_ylabel() or "").lower()
        title = (ax.get_title() or "").lower()
        is_trace_axis = "draw" in xlabel or "draw" in ylabel or "trace" in title
        is_distribution_axis = "ma_bmd" in title or "bmd" in xlabel or "density" in ylabel

        if is_trace_axis:
            for key in ("lower", "median", "upper"):
                ax.axhline(
                    stats[key],
                    color=colors[key],
                    linestyle=linestyles[key],
                    linewidth=1.25,
                    alpha=0.9,
                )
        elif is_distribution_axis:
            for key in ("lower", "median", "upper"):
                ax.axvline(
                    stats[key],
                    color=colors[key],
                    linestyle=linestyles[key],
                    linewidth=1.25,
                    alpha=0.9,
                    label=labels[key],
                )

    handles: list = []
    labels: list[str] = []
    for ax in fig.axes:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        handles.extend(ax_handles)
        labels.extend(ax_labels)

    if handles and fig.axes:
        by_label = dict(zip(labels, handles, strict=False))
        fig.axes[0].legend(by_label.values(), by_label.keys())


def _bmd_summary_table(idata: az.InferenceData, alpha: float) -> pd.DataFrame:
    lower_label = f"alpha_lower ({alpha:.2f})"
    upper_label = f"alpha_upper ({1 - alpha:.2f})"

    records: list[dict[str, float | str]] = []
    model_names = list(idata.posterior.coords["model"].values)
    bmd = np.asarray(idata.posterior["BMD"].values, dtype=float)

    for model_idx, model_name in enumerate(model_names):
        draws = bmd[:, :, model_idx].reshape(-1)
        records.append(
            {
                "model": str(model_name),
                "median": float(np.quantile(draws, 0.5)),
                lower_label: float(np.quantile(draws, alpha)),
                upper_label: float(np.quantile(draws, 1 - alpha)),
                "r_hat": np.nan,
                "ess_bulk": np.nan,
                "ess_tail": np.nan,
            }
        )

    ma_quantiles = _ma_bmd_quantiles(idata, alpha)
    records.append(
        {
            "model": "MA_BMD",
            "median": ma_quantiles["median"],
            lower_label: ma_quantiles["lower"],
            upper_label: ma_quantiles["upper"],
            "r_hat": np.nan,
            "ess_bulk": np.nan,
            "ess_tail": np.nan,
        }
    )

    return pd.DataFrame.from_records(records).set_index("model")


def _multi_summary_table(
    idata: az.InferenceData,
    var_names: list[str],
    hdi_prob: float,
) -> pd.DataFrame:
    def _fallback_summary(var_name: str) -> pd.DataFrame:
        values = np.asarray(idata.posterior[var_name].values, dtype=float)
        median_q = 0.5
        alpha = (1 - hdi_prob) / 2
        lower_q = alpha
        upper_q = 1 - alpha

        if values.ndim == 2:
            draws = values.reshape(-1)
            return pd.DataFrame(
                {
                    "median": [float(np.nanquantile(draws, median_q))],
                    f"eti_{lower_q:.0%}": [float(np.nanquantile(draws, lower_q))],
                    f"eti_{upper_q:.0%}": [float(np.nanquantile(draws, upper_q))],
                    "r_hat": [np.nan],
                    "ess_bulk": [np.nan],
                    "ess_tail": [np.nan],
                },
                index=pd.Index([var_name], name="var_name"),
            )

        if values.ndim == 3:
            rows = []
            labels = []
            model_labels = list(idata.posterior.coords["model"].values)
            for model_idx, model_name in enumerate(model_labels):
                draws = values[:, :, model_idx].reshape(-1)
                rows.append(
                    {
                        "median": float(np.nanquantile(draws, median_q)),
                        f"eti_{lower_q:.0%}": float(np.nanquantile(draws, lower_q)),
                        f"eti_{upper_q:.0%}": float(np.nanquantile(draws, upper_q)),
                        "r_hat": np.nan,
                        "ess_bulk": np.nan,
                        "ess_tail": np.nan,
                    }
                )
                labels.append(f"{var_name}[{model_idx}:{model_name}]")
            return pd.DataFrame(rows, index=pd.Index(labels, name="var_name"))

        raise ValueError(f"Unsupported posterior shape for {var_name}: {values.shape}")

    median_parts: list[pd.DataFrame] = []
    ess_parts: list[pd.DataFrame] = []

    for var_name in var_names:
        try:
            median_parts.append(
                az.summary(
                    idata,
                    var_names=[var_name],
                    hdi_prob=hdi_prob,
                    stat_focus="median",
                )
            )
            ess_parts.append(
                az.summary(
                    idata,
                    var_names=[var_name],
                    hdi_prob=hdi_prob,
                    stat_focus="mean",
                )
            )
        except ValueError:
            fallback = _fallback_summary(var_name)
            median_parts.append(fallback)
            ess_parts.append(fallback[["r_hat", "ess_bulk", "ess_tail"]])

    median_summary = pd.concat(median_parts, axis=0) if median_parts else pd.DataFrame()
    ess_summary = pd.concat(ess_parts, axis=0) if ess_parts else pd.DataFrame()

    median_summary = median_summary.drop(
        columns=[
            column for column in ("ess_median", "ess_tail", "ess_bulk") if column in median_summary
        ],
        errors="ignore",
    )
    ess_columns = [
        column for column in ("r_hat", "ess_bulk", "ess_tail") if column in ess_summary.columns
    ]
    if ess_columns:
        ess_to_join = ess_summary[ess_columns].drop(
            columns=[column for column in ess_columns if column in median_summary.columns],
            errors="ignore",
        )
        if len(ess_to_join.columns) > 0:
            median_summary = median_summary.join(ess_to_join)

    return median_summary


def _extract_model_row(summary: pd.DataFrame, var_name: str, model_name: str) -> pd.Series | None:
    for idx, row in summary.iterrows():
        label = str(idx)
        if label == var_name and model_name == "MA_BMD":
            return row
        if label == f"{var_name}[{model_name}]":
            return row
        if label.startswith(f"{var_name}[") and label.endswith(f":{model_name}]"):
            return row
    return None


def _drop_empty_summary_rows(df: pd.DataFrame) -> pd.DataFrame:
    # ArviZ can emit rows for parameter/model combinations that are entirely missing.
    return df.dropna(axis=0, how="all")


def _disttype_suffix(disttype: DistType | None) -> str:
    suffix_map = {
        DistType.normal: "CV",
        DistType.normal_ncv: "NCV",
        DistType.log_normal: "Lognormal",
    }
    return suffix_map.get(disttype, str(disttype))


def _model_color_map(model_names: list[str]) -> dict[str, tuple]:
    colors = plt.cm.tab10.colors
    return {model_name: colors[idx % len(colors)] for idx, model_name in enumerate(model_names)}


def _add_figure_legend(fig: plt.Figure, items: list[tuple[str, tuple]], ncol: int | None = None):
    handles = [Line2D([0], [0], color=color, lw=2, label=label) for label, color in items]
    fig.legend(
        handles=handles,
        labels=[label for label, _ in items],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=ncol or min(4, max(1, len(handles))),
        frameon=False,
    )


def _bmd_distributions_figure(idata: az.InferenceData) -> plt.Figure:
    model_names = [str(name) for name in idata.posterior.coords["model"].values]
    color_map = _model_color_map(model_names)
    fig, ax = plt.subplots(figsize=(8, 5))

    for model_name in model_names:
        draws = np.asarray(
            idata.posterior["BMD"].sel(model=model_name).values, dtype=float
        ).reshape(-1)
        draws = draws[np.isfinite(draws)]
        if draws.size == 0:
            continue
        az.plot_dist(draws, ax=ax, color=color_map[model_name])

    ma_draws = np.asarray(idata.posterior["MA_BMD"].values, dtype=float).reshape(-1)
    ma_draws = ma_draws[np.isfinite(ma_draws)]
    if ma_draws.size > 0:
        az.plot_dist(ma_draws, ax=ax, color="black", plot_kwargs={"linewidth": 2.5})

    ax.set_title("BMD distributions")
    ax.set_ylabel("Density")
    _add_figure_legend(
        fig,
        [(model_name, color_map[model_name]) for model_name in model_names]
        + [("Model Average", "black")],
        ncol=min(5, max(1, len(model_names) + 1)),
    )
    fig.tight_layout(rect=(0.03, 0.03, 0.97, 0.88))
    return fig


def _parameter_group_trace_figure(
    idata: az.InferenceData,
    group_name: str,
    model_names: list[str],
    param_names: list[str],
) -> plt.Figure:
    plot_names = ["BMD", *param_names]

    fig, axes = plt.subplots(
        len(plot_names),
        2,
        figsize=(11, 2.6 * len(plot_names)),
        squeeze=False,
    )
    color_map = _model_color_map(model_names)

    for row, var_name in enumerate(plot_names):
        ax_dist = axes[row, 0]
        ax_trace = axes[row, 1]

        for model_name in model_names:
            values = np.asarray(idata.posterior[var_name].sel(model=model_name).values, dtype=float)
            draws = values.reshape(-1)
            draws = draws[np.isfinite(draws)]
            if draws.size == 0:
                continue

            az.plot_dist(draws, ax=ax_dist, color=color_map[model_name])
            ax_trace.plot(draws, color=color_map[model_name], alpha=0.5, linewidth=0.8)

        ax_dist.set_title(var_name)
        ax_trace.set_title(var_name)

    _add_figure_legend(fig, [(model_name, color_map[model_name]) for model_name in model_names])
    fig.suptitle(f"{group_name} parameter distributions", y=0.985)
    fig.tight_layout(rect=(0.03, 0.03, 0.97, 0.965))
    return fig


def _parameter_group_records(
    idata: az.InferenceData,
    session,
    hdi_prob: float,
) -> list[dict[str, object]]:
    excluded_vars = {"BMD", "MA_BMD", "model_weights", "n_param"}

    grouped_models: dict[str, list] = {}
    for model in session.model_average.models:
        grouped_models.setdefault(model.bmd_model_class.verbose, []).append(model)

    records: list[dict[str, object]] = []
    for group_name, models in grouped_models.items():
        model_names = [model.name() for model in models]
        param_names: list[str] = []
        for model in models:
            for name in model.results.parameters.names:
                if name not in excluded_vars and name not in param_names:
                    param_names.append(name)

        rows: list[dict[str, object]] = []
        summary = _drop_empty_summary_rows(_multi_summary_table(idata, param_names, hdi_prob))
        for model_name in model_names:
            for param_name in param_names:
                stats = _extract_model_row(summary, param_name, model_name)
                if stats is None:
                    continue
                rows.append(
                    {
                        "Model": _disttype_suffix(
                            next(
                                model for model in models if model.name() == model_name
                            ).settings.disttype
                        ),
                        "Parameter": param_name,
                        **stats.to_dict(),
                    }
                )

        if not rows:
            continue

        figure = _parameter_group_trace_figure(idata, group_name, model_names, param_names)
        records.append(
            {
                "name": group_name,
                "model_names": model_names,
                "var_names": param_names,
                "summary": pd.DataFrame(rows),
                "trace_figure": figure,
            }
        )

    return records


def _bmd_diagnostics_table(idata: az.InferenceData, hdi_prob: float) -> pd.DataFrame:
    summary = _drop_empty_summary_rows(_multi_summary_table(idata, ["BMD", "MA_BMD"], hdi_prob))
    model_names = list(idata.posterior.coords["model"].values)
    rows: list[dict[str, object]] = []

    for model_name in model_names:
        stats = _extract_model_row(summary, "BMD", str(model_name))
        if stats is None:
            continue
        rows.append({"model": str(model_name), **stats.to_dict()})

    ma_stats = _extract_model_row(summary, "MA_BMD", "MA_BMD")
    if ma_stats is not None:
        rows.append({"model": "MA_BMD", **ma_stats.to_dict()})

    return pd.DataFrame(rows).set_index("model")


def get_model_average_figures(
    session,
    n_chains: int = 1,
) -> dict[str, plt.Figure | pd.DataFrame | az.InferenceData | float]:
    idata = model_average_to_inferencedata(session, n_chains=n_chains)
    out: dict[str, plt.Figure | pd.DataFrame | az.InferenceData | float] = {"idata": idata}

    alpha = session.model_average.models[0].settings.alpha
    out["alpha"] = alpha
    hdi_prob = 1 - (2 * alpha)
    out["hdi_prob"] = hdi_prob

    ma_bmd_quantiles = _ma_bmd_quantiles(idata, alpha)
    ma_bmd_hdi = _ma_bmd_hdi(idata, hdi_prob)
    out["ma_bmd_quantiles"] = ma_bmd_quantiles
    out["ma_bmd_hdi"] = ma_bmd_hdi

    bmd_summary = _bmd_diagnostics_table(idata, hdi_prob)
    out["bmd_summary"] = bmd_summary

    param_names = [
        v
        for v in idata.posterior.data_vars
        if v not in ["BMD", "MA_BMD", "model_weights", "n_param"]
    ]

    multi_var_names = ["BMD", "MA_BMD", *param_names]

    multi_summary = _drop_empty_summary_rows(_multi_summary_table(idata, multi_var_names, hdi_prob))
    out["multi_summary"] = multi_summary
    out["parameter_groups"] = _parameter_group_records(idata, session, hdi_prob)

    axes = az.plot_posterior(idata, var_names=["MA_BMD"])
    fig = _figure_from_axes(axes)
    _add_ma_bmd_reference_lines(
        fig,
        ma_bmd_hdi,
        lower_label="MA_BMD lower alpha HDI",
        upper_label="MA_BMD upper alpha HDI",
    )
    fig.tight_layout()
    out["posterior"] = fig
    out["overlay"] = _bmd_distributions_figure(idata)

    return out
