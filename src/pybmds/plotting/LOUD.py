from __future__ import annotations

from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from itertools import cycle


print("module:", az)
print("file:", getattr(az, "__file__", None))
print("version:", getattr(az, "__version__", None))
print("has InferenceData:", hasattr(az, "InferenceData"))
print("some names:", [x for x in dir(az) if "Inference" in x or "plot_" in x][:20])


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


def model_average_to_inferencedata(session, n_chains: int = 1) -> az.InferenceData:
    """
    Convert a pybmds session with continuous model averaging results into an
    ArviZ InferenceData object.

    Posterior group:
        - one variable per pybmds parameter name (for example m0, m1, n, Var0)
        - BMD           : model-specific BMD draws
        - MA_BMD        : model-averaged BMD draws
        - model_weights : posterior model weights
        - n_param       : number of parameters per model

    Observed data group:
        - built from session.dataset
    """
    ma_result = session.model_average.results
    models = session.models

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

    return az.InferenceData(
        posterior=model_coord_posterior,
        observed_data=observed_data,
    )


def save_model_average_inferencedata(
    session,
    path: str | Path,
    n_chains: int = 1,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    idata = model_average_to_inferencedata(session, n_chains=n_chains)
    idata.to_netcdf(path)
    return path


def get_model_average_figures(
    session,
    n_chains: int = 1,
) -> dict[str, plt.Figure]:
    idata = model_average_to_inferencedata(session, n_chains=n_chains)
    out: dict[str, plt.Figure] = {}

    # Posterior plot
    axes = az.plot_posterior(idata, var_names=["MA_BMD"])
    fig = np.ravel(np.atleast_1d(axes))[0].figure
    out["posterior"] = fig

    # Trace plot
    axes = az.plot_trace(idata, var_names=["MA_BMD"])
    fig = np.ravel(np.atleast_1d(axes))[0].figure
    out["trace"] = fig

    # Overlay distribution plot
    fig, ax = plt.subplots(figsize=(8, 5))

    # create a color cycle
    colors = plt.cm.tab10.colors  # nice distinct colors
    color_cycle = cycle(colors)

    for model_name in idata.posterior.coords["model"].values:
        color = next(color_cycle)
        az.plot_dist(
            idata.posterior["BMD"].sel(model=model_name),
            ax=ax,
            label=str(model_name),
            color=color,
        )

    # model average in bold / standout color
    az.plot_dist(
        idata.posterior["MA_BMD"],
        ax=ax,
        label="Model Average",
        color="black",
        plot_kwargs={"lw": 3},
    )

    ax.legend()
    ax.set_title("BMD distributions")
    out["overlay"] = fig

    return out
