import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy.special import gamma
from scipy.stats import norm

from ..session import Session
from .common import create_empty_figure


def _calc_alpha(row, dose_to_phi) -> float | None:
    """
    For ilc+ models, if phi is estimated for a specific dose group, we calculate the parameters of a beta distribution.
    This function calculates alpha.
    """
    if row["ilc_est"] == 1:
        phi_col = dose_to_phi.get(row["dose"])
        if phi_col and phi_col in row and row[phi_col] != 0:
            return row["est_prob"] * (1 / row[phi_col] - 1)
    return None


def _calc_rd(row, dose_to_phi) -> float | None:
    """
    For ilc+ models with phi_i > 0, we calculate Rd. For x in (0, 1, 2, …, n) the probability mass function
    defining the probability of observing x responders in a litter of size n in group d is defined as Rd.
    """
    if row["ilc_est"] == 1:
        alpha = row["alpha"]
        beta = row["beta"]
        litter_size = row["litter_size"]

        phi_column = dose_to_phi.get(row["dose"])

        if phi_column and phi_column in row and row[phi_column] == 0:
            return None

        if any(x is None or x <= 0 or math.isnan(x) for x in [alpha, beta, litter_size]):
            return None

        return (
            gamma(litter_size + 1)
            * gamma(alpha + beta)
            / (gamma(alpha) * gamma(beta) * gamma(litter_size + alpha + beta))
        )

    return None


def _generate_litter_columns(df: pd.DataFrame, dose_to_phi: dict) -> pd.DataFrame:
    """
    This function goes through the methodology for both ilc+ and ilc- models.

    If an ilc+ model has phi = 0, it will follow the ilc- approach.

    * "litter" is the binomial coefficient representing the number of ways to select (i)
        responders from a litter of size (n). In ilc+ models, this term is always 1,
        as the litter-level variability is already incorporated into the beta-binomial
        structure. In ilc- models, this follows the standard binomial formulation.
    * "litter_next" is the conditional probability of exactly (i) responders in a litter
        of size (n). In ilc+ models, this is computed using the beta-binomial probability.
        In ilc- models, this is computed using the binomial probability.
    * "litter_final" is the expected number of litters with exactly (i) responders in dose
        group (d), obtained by multiplying litter and litter_next.
        For ilc+ models, since litter=1, we have litter_final=litter_next*Rd.
        For ilc- models, we have litter_final=litter*litter_next.
    """
    max_litter_size = int(df["litter_size"].max())

    litter = pd.DataFrame(index=df.index)
    litter_next = pd.DataFrame(index=df.index)
    litter_final = pd.DataFrame(index=df.index)

    for i in range(max_litter_size + 1):

        def calc_litter(row, i=i):
            phi_column = dose_to_phi.get(row["dose"], None)
            is_ilc_plus = row["ilc_est"] == 1
            phi_is_zero = phi_column and phi_column in row and row[phi_column] == 0

            if is_ilc_plus and not phi_is_zero:
                return 1
            else:
                if i > int(row["litter_size"]):
                    return 0
                return math.comb(int(row["litter_size"]), i)

        def calc_litter_next(row, i=i):
            phi_column = dose_to_phi.get(row["dose"], None)
            is_ilc_plus = row["ilc_est"] == 1
            phi_is_zero = phi_column and phi_column in row and row[phi_column] == 0

            if is_ilc_plus and not phi_is_zero:
                if i > row["litter_size"]:
                    return 0
                return (
                    gamma(i + row["alpha"])
                    * gamma(row["litter_size"] - i + row["beta"])
                    / (gamma(i + 1) * gamma(row["litter_size"] - i + 1))
                )
            else:
                return row["est_prob"] ** i * (1 - row["est_prob"]) ** (row["litter_size"] - i)

        def calc_litter_final(row, i=i):
            phi_column = dose_to_phi.get(row["dose"], None)
            is_ilc_plus = row["ilc_est"] == 1
            phi_is_zero = phi_column and phi_column in row and row[phi_column] == 0

            base = litter.at[row.name, i] * litter_next.at[row.name, i]
            return base * row["Rd"] if is_ilc_plus and not phi_is_zero else base

        litter[i] = df.apply(calc_litter, axis=1)
        litter_next[i] = df.apply(calc_litter_next, axis=1)
        litter_final[i] = df.apply(calc_litter_final, axis=1)

    return pd.concat([df.copy(), litter_final], axis=1)


def _sum_litter_by_dose(df: pd.DataFrame) -> pd.DataFrame:
    """
    Estimated counts for each dose group.
    """
    litter_final_cols = [col for col in df.columns if isinstance(col, int)]
    return df.groupby(["dose", "model_name"])[litter_final_cols].sum().reset_index()


def _count_litter_obs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the number of litters in each dose group that have exactly x responders
    """
    litter_final_cols = [col for col in df.columns if isinstance(col, int)]

    grouped_df = df.groupby(["dose", "model_name"])
    result_df = pd.DataFrame(
        data=dict(
            dose=grouped_df["dose"].first().values,
            model_name=grouped_df["model_name"].first().values,
        )
    )
    for litter_col in litter_final_cols:
        counts = []
        for _, group in grouped_df:
            count = (group["obs"] == litter_col).sum()
            counts.append(count)
        result_df[litter_col] = counts
    result_df.reset_index(drop=True, inplace=True)
    return result_df


def _add_all_doses(counted_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds rows at the end for the estimated and observed values for all doses
    """
    summed_df = (
        counted_df.groupby("model_name")[counted_df.columns.difference(["dose", "model_name"])]
        .sum()
        .reset_index()
        .assign(dose="All")
    )
    return pd.concat([counted_df, summed_df], ignore_index=True)


def _transform(observed_df: pd.DataFrame, estimated_df: pd.DataFrame) -> pd.DataFrame:
    obs_long = observed_df.melt(
        id_vars=["dose", "model_name"], var_name="litter_final", value_name="obs_count"
    )
    est_long = estimated_df.melt(
        id_vars=["dose", "model_name"], var_name="litter_final", value_name="est_count"
    )
    return pd.merge(obs_long, est_long, on=["dose", "model_name", "litter_final"], how="outer")


def _calculate_plotting(incidence: float, n: float, alpha: float = 0.05):
    """
    Add confidence intervals to dichotomous datasets.
    https://www.epa.gov/sites/production/files/2020-09/documents/bmds_3.2_user_guide.pdf

    The error bars shown in BMDS plots use alpha = 0.05 and so represent the 95% confidence
    intervals on the observed proportions (independent of model).
    """
    p = incidence / float(n)
    z = norm.ppf(1 - alpha / 2)
    z2 = z * z
    q = 1.0 - p
    tmp1 = 2 * n * p + z2
    ll = ((tmp1 - 1) - z * np.sqrt(z2 - (2 + 1 / n) + 4 * p * (n * q + 1))) / (2 * (n + z2))
    ul = ((tmp1 + 1) + z * np.sqrt(z2 + (2 + 1 / n) + 4 * p * (n * q - 1))) / (2 * (n + z2))
    return p, ll, ul


def dose_litter_response_data(session: Session) -> pd.DataFrame:
    """
    Generate data for plotting nested litter plots comparing observed vs estimated responder
    distributions for each dose level.

    We extract the necessary results and parameter values to create the plots and save
    them in "final". We need the model parameters to calculate alpha and beta and the
    "Litter Data" information such as the estimated probabilities, litter size, etc. for
    each litter. The model parameters and phi's are saved as columns and each value repeated
    for each bmds_model_index for individual litter calculations.
    """
    model_res = []
    for j, model in enumerate(session.models):
        if not model.has_results:
            continue

        res = model.results
        name = model.name()
        param_map = dict(zip(res.parameter_names, res.parameters, strict=True))
        ilc_est = model.settings.intralitter_correlation
        for i in range(len(res.litter.dose)):
            model_res.append(
                {
                    "bmds_model_index": j,
                    "model_name": name,
                    "ilc_est": ilc_est,
                    "dose": res.litter.dose[i],
                    "est_prob": res.litter.estimated_probabilities[i],
                    "lsc": res.litter.lsc[i],
                    "litter_size": res.litter.litter_size[i],
                    "obs": res.litter.observed[i],
                    "scaled_residuals": res.litter.scaled_residuals[i],
                    **param_map,
                }
            )

    final = pd.DataFrame(model_res)

    unique_doses = sorted(final["dose"].unique())
    dose_to_phi = {dose: f"phi{index + 1}" for index, dose in enumerate(unique_doses)}

    final["alpha"] = final.apply(lambda row: _calc_alpha(row, dose_to_phi), axis=1)
    final["beta"] = final["alpha"] * (1 - final["est_prob"]) / final["est_prob"]
    final["Rd"] = final.apply(lambda row: _calc_rd(row, dose_to_phi), axis=1)

    df = _generate_litter_columns(final, dose_to_phi)
    estimated = _add_all_doses(_sum_litter_by_dose(df))
    observed = _add_all_doses(_count_litter_obs(df))
    plot_data = _transform(observed, estimated)

    plot_data["total_litters_per_dose"] = plot_data.groupby(["dose", "model_name"])[
        "obs_count"
    ].transform("sum")
    plot_data[["obs_proportion", "obs_ci_lower", "obs_ci_upper"]] = plot_data.apply(
        lambda row: _calculate_plotting(row["obs_count"], row["total_litters_per_dose"]),
        axis=1,
        result_type="expand",
    )
    plot_data["obs_ci_lower_count"] = plot_data.obs_ci_lower * plot_data.total_litters_per_dose
    plot_data["obs_ci_upper_count"] = plot_data.obs_ci_upper * plot_data.total_litters_per_dose
    return plot_data


def dose_litter_response_plot(session: Session) -> Figure:
    """
    Generate nested litter plots comparing observed vs estimated responder distributions
    for each dose level.

    Args:
        session (pybmds.Session): A fitted BMDS session with nested models.

    Returns:
        Figure: a maptlotlib figure.
    """
    plot_data = dose_litter_response_data(session)
    plot_data["lls"] = (plot_data.obs_count - plot_data.obs_ci_lower_count).clip(lower=0)
    plot_data["uls"] = (plot_data.obs_ci_upper_count - plot_data.obs_count).clip(lower=0)

    unique_doses = plot_data.dose.unique()
    models = plot_data.model_name.unique()
    nrows = plot_data.dose.nunique()

    height = 2 * nrows
    width = plot_data.litter_final.max() / 2

    fig = create_empty_figure(rows=nrows, cols=1, figsize=(width, height))
    cmap = plt.get_cmap("tab10")
    model_colors = {model: cmap(i % cmap.N) for i, model in enumerate(models)}

    for ax, dose in zip(fig.get_axes(), unique_doses, strict=True):
        dose_data = plot_data[plot_data.dose == dose]
        for model in models:
            model_data = dose_data.query(f"model_name == '{model}'")
            color = model_colors[model]
            ax.errorbar(
                model_data.litter_final,
                model_data.obs_count,
                yerr=(model_data.lls, model_data.uls),
                fmt="o",
                color="grey",
                capsize=4,
                markersize=5,
            )
            ax.plot(
                model_data.litter_final,
                model_data.est_count,
                color=color,
                marker=".",
                linestyle="-",
            )
        ax.margins(x=0.05, y=0.05)
        ax.set_title(f"Dose = {dose}")
        ax.set_xlabel("# of Responders")
        ax.set_ylabel("# of Litters")
        ax.tick_params(axis="x")
        ax.set_xticks(range(int(plot_data["litter_final"].max()) + 1))

    # build legend
    legend = [
        Line2D([0], [0], marker="o", color=model_colors[model_name], label=f"{model_name}")
        for model_name in models
    ]
    obs = Line2D([0], [0], marker="o", color="grey", linestyle="None", label="Observations")
    legend.insert(0, obs)

    fig.legend(
        legend,
        [h.get_label() for h in legend],
        loc="upper right",
        bbox_to_anchor=(1.0, 0.96),
    )

    # leave space for legend
    fig.tight_layout(rect=[0, 0, 0.7, 1])
    return fig
