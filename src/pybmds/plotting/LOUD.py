from __future__ import annotations

from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from ..constants import DistType

az._log.disabled = True


def _reshape_draws(arr: np.ndarray, n_chains: int) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = arr.shape[0]
    if n % n_chains != 0:
        raise ValueError(f"{n=} must be divisible by {n_chains=}.")
    draws_per_chain = n // n_chains
    return arr.reshape((n_chains, draws_per_chain) + arr.shape[1:])


def _nan_nonfinite(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=float).copy()
    arr[~np.isfinite(arr)] = np.nan
    return arr


def _pad_draws(arr: np.ndarray, n_draws: int) -> np.ndarray:
    arr = _nan_nonfinite(arr)
    if arr.ndim == 1:
        out = np.full(n_draws, np.nan, dtype=float)
        rows = min(arr.shape[0], n_draws)
        out[:rows] = arr[:rows]
        return out

    if arr.ndim == 2:
        out = np.full((n_draws, arr.shape[1]), np.nan, dtype=float)
        rows = min(arr.shape[0], n_draws)
        out[:rows, :] = arr[:rows, :]
        return out

    raise ValueError(f"Unsupported draw shape: {arr.shape}")


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

    ma_bmd = _nan_nonfinite(ma_result.bmd_dist)
    model_bmd_draws = [_nan_nonfinite(bmds) for bmds in ma_result.model_bmd_dist]
    model_parm_draws = [_nan_nonfinite(parms) for parms in ma_result.model_parm_dist]

    n_draws = max([len(ma_bmd)] + [len(bmds) for bmds in model_bmd_draws])
    n_models = len(model_names)
    n_params = len(all_param_names)

    model_bmd = np.full((n_draws, n_models), np.nan, dtype=float)
    model_params = np.full((n_draws, n_models, n_params), np.nan, dtype=float)
    n_param_by_model = np.zeros(n_models, dtype=int)

    for model_idx, (bmds, parms, pnames) in enumerate(
        zip(model_bmd_draws, model_parm_draws, param_names_by_model, strict=True)
    ):
        if parms.ndim != 2:
            raise ValueError(
                f"Model {model_names[model_idx]!r} parameter draws must be 2D; got shape {parms.shape}."
            )

        if len(bmds) != parms.shape[0]:
            raise ValueError(
                f"Model {model_names[model_idx]!r} has {len(bmds)} BMD draws but "
                f"{parms.shape[0]} parameter draws."
            )

        if parms.shape[1] != len(pnames):
            raise ValueError(
                f"Model {model_names[model_idx]!r} has {parms.shape[1]} parameters but "
                f"{len(pnames)} names."
            )

        bmds = _pad_draws(bmds, n_draws)
        parms = _pad_draws(parms, n_draws)
        model_bmd[:, model_idx] = bmds
        n_param_by_model[model_idx] = len(pnames)

        for param_idx, pname in enumerate(pnames):
            global_idx = all_param_names.index(pname)
            model_params[:, model_idx, global_idx] = parms[:, param_idx]

    ma_bmd = _pad_draws(ma_bmd, n_draws)
    posteriors = _nan_nonfinite(ma_result.posteriors)

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
                _reshape_draws(ma_bmd, n_chains),
            ),
            "weights": (
                ("model",),
                posteriors,
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
    ma_bmd = ma_bmd[np.isfinite(ma_bmd)]
    return {
        "lower": float(np.nanquantile(ma_bmd, alpha)),
        "median": float(np.nanquantile(ma_bmd, 0.5)),
        "upper": float(np.nanquantile(ma_bmd, 1 - alpha)),
    }


def _ma_bmd_hdi(idata: az.InferenceData, hdi_prob: float) -> dict[str, float]:
    ma_bmd = np.asarray(idata.posterior["MA_BMD"].values, dtype=float).reshape(-1)
    ma_bmd = ma_bmd[np.isfinite(ma_bmd)]
    lower, upper = np.asarray(az.hdi(ma_bmd, hdi_prob=hdi_prob), dtype=float)
    return {
        "lower": float(lower),
        "median": float(np.nanquantile(ma_bmd, 0.5)),
        "upper": float(upper),
    }


def _ma_bmd_posterior_figure(
    idata: az.InferenceData,
    stats: dict[str, float],
) -> plt.Figure:
    ma_draws = np.asarray(idata.posterior["MA_BMD"].values, dtype=float).reshape(-1)
    ma_draws = ma_draws[np.isfinite(ma_draws)]

    fig, ax = plt.subplots(figsize=(8, 5))
    if ma_draws.size > 0:
        az.plot_dist(ma_draws, ax=ax, color="black", plot_kwargs={"linewidth": 2.5})

    ax.axvline(
        stats["lower"],
        color="#B04A3A",
        linestyle="--",
        linewidth=1.5,
        alpha=0.9,
        label="BMDL",
    )
    ax.axvline(
        stats["upper"],
        color="#2A7F62",
        linestyle="--",
        linewidth=1.5,
        alpha=0.9,
        label="BMDU",
    )
    ax.axvline(
        stats["median"],
        color="#111111",
        linestyle="-",
        linewidth=1.5,
        alpha=0.9,
        label="BMD",
    )
    text_y_positions = {"lower": 0.82, "median": 0.90, "upper": 0.74}
    labels = {"lower": "BMDL", "median": "BMD", "upper": "BMDU"}
    for key in ("lower", "median", "upper"):
        ax.annotate(
            f"{labels[key]}: {stats[key]:.4g}",
            xy=(stats[key], text_y_positions[key]),
            xycoords=("data", "axes fraction"),
            xytext=(6, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=10,
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "alpha": 0.85,
                "edgecolor": "none",
            },
        )
    ax.set_title("Model-averaged BMD distribution", pad=6)
    ax.set_xlabel("BMD")
    ax.set_ylabel("Density")
    ax.legend()
    fig.tight_layout()
    return fig


def _bmd_summary_table(idata: az.InferenceData, alpha: float) -> pd.DataFrame:
    records: list[dict[str, float | str]] = []
    model_names = list(idata.posterior.coords["model"].values)
    bmd = np.asarray(idata.posterior["BMD"].values, dtype=float)

    for model_idx, model_name in enumerate(model_names):
        draws = bmd[:, :, model_idx].reshape(-1)
        draws = draws[np.isfinite(draws)]
        records.append(
            {
                "model": str(model_name),
                "BMD": float(np.nanquantile(draws, 0.5)),
                "BMDL": float(np.nanquantile(draws, alpha)),
                "BMDU": float(np.nanquantile(draws, 1 - alpha)),
                "r_hat": np.nan,
                "ess_bulk": np.nan,
                "ess_tail": np.nan,
            }
        )

    ma_quantiles = _ma_bmd_quantiles(idata, alpha)
    records.append(
        {
            "model": "MA_BMD",
            "BMD": ma_quantiles["median"],
            "BMDL": ma_quantiles["lower"],
            "BMDU": ma_quantiles["upper"],
            "r_hat": np.nan,
            "ess_bulk": np.nan,
            "ess_tail": np.nan,
        }
    )

    return _rename_summary_columns(pd.DataFrame.from_records(records).set_index("model"))


def _percent_label_from_eti_column(column: str) -> str:
    return column.removeprefix("eti_")


def _rename_summary_columns(summary: pd.DataFrame, bmd_labels: bool = False) -> pd.DataFrame:
    rename: dict[str, str] = {}
    if "median" in summary.columns:
        rename["median"] = "BMD" if bmd_labels else "Median"

    readable_labels = {
        "mad": "Median Absolute Deviation",
        "mcse_median": "Markov Chain Standard Error (Median)",
        "r_hat": "R-hat",
        "ess_bulk": "Bulk Effective Sample Size",
        "ess_tail": "Tail Effective Sample Size",
    }
    rename.update(
        {column: label for column, label in readable_labels.items() if column in summary.columns}
    )

    eti_columns = [column for column in summary.columns if column.startswith("eti_")]
    if len(eti_columns) >= 2:
        eti_columns = sorted(
            eti_columns,
            key=lambda column: float(column.removeprefix("eti_").removesuffix("%")),
        )
        if bmd_labels:
            rename[eti_columns[0]] = "BMDL"
            rename[eti_columns[-1]] = "BMDU"
        else:
            rename.update(
                {column: _percent_label_from_eti_column(column) for column in eti_columns}
            )

    return summary.rename(columns=rename)


def _summary_from_draws(
    draws: np.ndarray, var_name: str, label: str, hdi_prob: float
) -> pd.DataFrame:
    draws = np.asarray(draws, dtype=float).reshape(-1)
    draws = draws[np.isfinite(draws)]
    if draws.size == 0:
        return pd.DataFrame()

    dataset = xr.Dataset(
        {var_name: (("chain", "draw"), draws.reshape(1, draws.size))},
        coords={"chain": [0], "draw": np.arange(draws.size)},
    )
    idata = az.InferenceData(posterior=dataset)

    try:
        median_summary = az.summary(
            idata,
            var_names=[var_name],
            hdi_prob=hdi_prob,
            stat_focus="median",
        )
        ess_summary = az.summary(
            idata,
            var_names=[var_name],
            hdi_prob=hdi_prob,
            stat_focus="mean",
        )
        median_summary = median_summary.iloc[[0]].copy()
        ess_summary = ess_summary.iloc[[0]].copy()
        median_summary.index = pd.Index([label], name="var_name")
        ess_summary.index = pd.Index([label], name="var_name")
        median_summary = median_summary.drop(
            columns=[
                column
                for column in ("ess_median", "ess_tail", "ess_bulk")
                if column in median_summary
            ],
            errors="ignore",
        )
        ess_columns = [
            column for column in ("r_hat", "ess_bulk", "ess_tail") if column in ess_summary.columns
        ]
        if ess_columns:
            median_summary = median_summary.join(
                ess_summary[ess_columns].drop(
                    columns=[column for column in ess_columns if column in median_summary.columns],
                    errors="ignore",
                )
            )
        return median_summary
    except ValueError:
        alpha = (1 - hdi_prob) / 2
        return pd.DataFrame(
            {
                "median": [float(np.nanquantile(draws, 0.5))],
                f"eti_{alpha:.0%}": [float(np.nanquantile(draws, alpha))],
                f"eti_{1 - alpha:.0%}": [float(np.nanquantile(draws, 1 - alpha))],
                "r_hat": [np.nan],
                "ess_bulk": [np.nan],
                "ess_tail": [np.nan],
            },
            index=pd.Index([label], name="var_name"),
        )


def _finite_summary_table(idata: az.InferenceData, var_names: list[str], hdi_prob: float):
    parts: list[pd.DataFrame] = []
    model_labels = (
        list(idata.posterior.coords["model"].values) if "model" in idata.posterior.coords else []
    )

    for var_name in var_names:
        values = np.asarray(idata.posterior[var_name].values, dtype=float)
        if values.ndim == 2:
            parts.append(_summary_from_draws(values, var_name, var_name, hdi_prob))
        elif values.ndim == 3:
            for model_idx, model_name in enumerate(model_labels):
                label = f"{var_name}[{model_idx}:{model_name}]"
                parts.append(
                    _summary_from_draws(values[:, :, model_idx], var_name, label, hdi_prob)
                )
        else:
            raise ValueError(f"Unsupported posterior shape for {var_name}: {values.shape}")

    parts = [part for part in parts if not part.empty]
    return pd.concat(parts, axis=0) if parts else pd.DataFrame()


def _multi_summary_table(
    idata: az.InferenceData,
    var_names: list[str],
    hdi_prob: float,
) -> pd.DataFrame:
    if idata is not None:
        return _finite_summary_table(idata, var_names, hdi_prob)

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


def _parameter_group_model_label(model) -> str:
    disttype = getattr(model.settings, "disttype", None)
    if disttype is not None:
        return _disttype_suffix(disttype)
    return model.name()


def _model_color_map(model_names: list[str]) -> dict[str, tuple]:
    bmds_colors = {
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
    colors = plt.cm.tab10.colors
    return {
        model_name: bmds_colors.get(model_name, colors[idx % len(colors)])
        for idx, model_name in enumerate(model_names)
    }


def _add_figure_legend(fig: plt.Figure, items: list[tuple[str, tuple]], ncol: int | None = None):
    handles = [Line2D([0], [0], color=color, lw=2, label=label) for label, color in items]
    fig.legend(
        handles=handles,
        labels=[label for label, _ in items],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=ncol or min(3, max(1, len(handles))),
        frameon=False,
    )


def _bmd_distributions_figure(idata: az.InferenceData) -> plt.Figure:
    model_names = [str(name) for name in idata.posterior.coords["model"].values]
    color_map = _model_color_map(model_names)
    fig, ax = plt.subplots(figsize=(8, 6.5))

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

    ax.set_title("BMD distributions", pad=6)
    ax.set_ylabel("Density")
    _add_figure_legend(
        fig,
        [(model_name, color_map[model_name]) for model_name in model_names]
        + [("Model Average", "black")],
        ncol=min(5, max(1, len(model_names) + 1)),
    )
    fig.tight_layout(rect=(0.03, 0.03, 0.97, 0.80))
    return fig


def _parameter_figure_height(n_rows: int) -> float:
    row_height = 2.5 if n_rows <= 4 else 2.25 if n_rows <= 7 else 2
    return min(11.0, max(8.0, row_height * n_rows))


def _parameter_group_trace_figure(
    idata: az.InferenceData,
    group_name: str,
    model_names: list[str],
    param_names: list[str],
    color_map: dict[str, tuple | str],
    param_model_names: dict[str, set[str]] | None = None,
) -> plt.Figure:
    plot_names = ["BMD", *param_names]

    fig = Figure(figsize=(11, _parameter_figure_height(len(plot_names))))
    axes = fig.subplots(len(plot_names), 2, squeeze=False)
    for row, var_name in enumerate(plot_names):
        ax_dist = axes[row, 0]
        ax_trace = axes[row, 1]

        for model_name in model_names:
            if var_name != "BMD" and param_model_names is not None:
                if model_name not in param_model_names.get(var_name, set()):
                    continue
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
    compressed: bool = True,
) -> list[dict[str, object]]:
    excluded_vars = {"BMD", "MA_BMD", "model_weights", "n_param"}
    all_model_names = [model.name() for model in session.model_average.models]
    color_map = _model_color_map(all_model_names)

    grouped_models: dict[str, list] = {}
    for model in session.model_average.models:
        if compressed:
            group_name = model.bmd_model_class.verbose
        else:
            group_name = f"{model.bmd_model_class.verbose} {_parameter_group_model_label(model)}"
        grouped_models.setdefault(group_name, []).append(model)

    records: list[dict[str, object]] = []
    for group_name, models in grouped_models.items():

        def add_record(record_name: str, record_models: list, include_param) -> None:
            model_names = [model.name() for model in record_models]
            param_names: list[str] = []
            param_model_names: dict[str, set[str]] = {}
            for model in record_models:
                model_name = model.name()
                for name in model.results.parameters.names:
                    if name in excluded_vars or not include_param(model, name):
                        continue
                    if name not in param_names:
                        param_names.append(name)
                    param_model_names.setdefault(name, set()).add(model_name)

            if not param_names:
                return

            rows: list[dict[str, object]] = []
            summary = _drop_empty_summary_rows(_multi_summary_table(idata, param_names, hdi_prob))
            summary = _rename_summary_columns(summary)
            for model in record_models:
                model_name = model.name()
                for param_name in param_names:
                    if model_name not in param_model_names.get(param_name, set()):
                        continue
                    stats = _extract_model_row(summary, param_name, model_name)
                    if stats is None:
                        continue
                    rows.append(
                        {
                            "Model": _parameter_group_model_label(model),
                            "Parameter": param_name,
                            **stats.to_dict(),
                        }
                    )

            if not rows:
                return

            figure = _parameter_group_trace_figure(
                idata, record_name, model_names, param_names, color_map, param_model_names
            )
            records.append(
                {
                    "name": record_name,
                    "model_names": model_names,
                    "var_names": param_names,
                    "summary": pd.DataFrame(rows),
                    "trace_figure": figure,
                }
            )

        add_record(
            group_name,
            models,
            lambda model, name: True,
        )

    return records


def _bmd_diagnostics_table(idata: az.InferenceData, hdi_prob: float) -> pd.DataFrame:
    summary = _drop_empty_summary_rows(_multi_summary_table(idata, ["BMD", "MA_BMD"], hdi_prob))
    summary = _rename_summary_columns(summary, bmd_labels=True)
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
    compressed: bool = True,
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
    out["parameter_groups"] = _parameter_group_records(
        idata, session, hdi_prob, compressed=compressed
    )

    out["posterior"] = _ma_bmd_posterior_figure(idata, ma_bmd_quantiles)
    out["overlay"] = _bmd_distributions_figure(idata)

    return out
