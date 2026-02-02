"""
Rao-Scott Transformation, original citation:

Fox JF, Hogan KA, Davis A. Dose-Response Modeling with Summary Data from Developmental Toxicity
Studies. Risk Anal. 2017 May;37(5):905-917. doi: 10.1111/risa.12667. PMID: 27567129.
"""

from enum import StrEnum
from io import BytesIO
from typing import ClassVar, NamedTuple

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from ... import plotting
from ...reporting.styling import Report, add_mpl_figure, df_to_table, write_setting_p
from ..dichotomous import DichotomousDataset


class Species(StrEnum):
    rat = "rat"
    mouse = "mouse"
    rabbit = "rabbit"


class Regression(StrEnum):
    least_square = "Least Squares"
    orthogonal = "Orthogonal"


class P(NamedTuple):
    a: float
    b: float
    sigma: float


class RaoScott:
    # see Table 2 from paper
    # fmt: off
    adjustment_parameters: ClassVar = {
        (Species.mouse,  Regression.least_square): P(1.5938, 0.2866, 0.2078),
        (Species.mouse,  Regression.orthogonal):   P(1.6943, 0.3132, 0.1863),
        (Species.rat,    Regression.least_square): P(1.6852, 0.3310, 0.1248),
        (Species.rat,    Regression.orthogonal):   P(1.8327, 0.3690, 0.1090),
        (Species.rabbit, Regression.least_square): P(1.0582, 0.2397, 0.1452),
        (Species.rabbit, Regression.orthogonal):   P(1.1477, 0.2739, 0.1299),
    }
    # fmt: on

    def __init__(self, dataset: DichotomousDataset, species: Species):
        self.dataset = dataset
        self.species = species
        self.df = self.calculate_design_effects_with_limit()

    def calculate_design_effects_with_limit(self) -> pd.DataFrame:
        df = pd.DataFrame(
            {
                "dose": self.dataset.doses,
                "incidence": self.dataset.incidences,
                "n": self.dataset.ns,
            }
        )
        df["fraction_affected"] = df.incidence / df.n

        # Masks for zero and nonzero fraction_affected
        zero_mask = df["fraction_affected"] == 0
        nonzero_mask = ~zero_mask

        # Prepare columns with NaN
        df["design_ls"] = np.nan
        df["design_o"] = np.nan
        df["design_avg"] = np.nan

        # Calculate design effects for nonzero fraction_affected
        if nonzero_mask.any():
            p_ls = self.adjustment_parameters[(self.species, Regression.least_square)]
            p_or = self.adjustment_parameters[(self.species, Regression.orthogonal)]

            log_fa = np.log(df.loc[nonzero_mask, "fraction_affected"])
            df.loc[nonzero_mask, "design_ls"] = np.exp(
                p_ls.a + (p_ls.b * log_fa) + (0.5 * p_ls.sigma)
            )
            df.loc[nonzero_mask, "design_o"] = np.exp(
                p_or.a + (p_or.b * log_fa) + (0.5 * p_or.sigma)
            )
            df.loc[nonzero_mask, "design_avg"] = df.loc[
                nonzero_mask, ["design_ls", "design_o"]
            ].mean(axis=1)

        # For zero_mask, set design_ls, design_o, and design_avg to 1
        df.loc[zero_mask, "design_ls"] = 1
        df.loc[zero_mask, "design_o"] = 1
        df.loc[zero_mask, "design_avg"] = 1

        # Detect and identify dose groups where the original design effect was less than 1
        df["design_ls_below_1"] = (df["design_ls"] < 1) & df["design_ls"].notna()
        df["design_o_below_1"] = (df["design_o"] < 1) & df["design_o"].notna()
        df["design_avg_below_1"] = (df["design_avg"] < 1) & df["design_avg"].notna()

        # Apply lower limit of 1 to the design effects
        df["design_ls"] = df["design_ls"].clip(lower=1)
        df["design_o"] = df["design_o"].clip(lower=1)
        df["design_avg"] = df["design_avg"].clip(lower=1)

        # Calculate adjusted columns using design_avg
        df["incidence_adjusted"] = df["incidence"] / df["design_avg"]
        df["n_adjusted"] = df["n"] / df["design_avg"]

        return df

    def figure(self, figsize: tuple[float, float] | None = None) -> Figure:
        fig = plotting.create_empty_figure(rows=1, cols=2, figsize=figsize)
        ax1, ax2 = fig.axes

        # N
        ax1.set_title("Original N vs Adjusted N")
        ax1.set_xlabel("Dose")
        ax1.set_ylabel("N")
        ax1.margins(plotting.PLOT_MARGINS)
        ax1.plot(
            self.df.dose,
            self.df.n,
            "o-",
            color="FireBrick",
            label="Original N",
            markersize=8,
            markeredgewidth=1,
            markeredgecolor="white",
        )
        ax1.plot(
            self.df.dose,
            self.df.n_adjusted,
            "^-",
            color="LightSalmon",
            label="Adjusted N",
            markersize=8,
            markeredgewidth=1,
            markeredgecolor="white",
        )
        legend = ax1.legend(**plotting.LEGEND_OPTS)
        for handle in legend.legend_handles:
            handle.set_markersize(8)

        # Incidence
        ax2.set_title("Original Incidence vs Adjusted Incidence")
        ax2.set_xlabel("Dose")
        ax2.set_ylabel("Incidence")
        ax2.margins(plotting.PLOT_MARGINS)
        ax2.plot(
            self.df.dose,
            self.df.incidence,
            "o-",
            color="MidnightBlue",
            label="Original Incidence",
            markersize=8,
            markeredgewidth=1,
            markeredgecolor="white",
        )
        ax2.plot(
            self.df.dose,
            self.df.incidence_adjusted,
            "^-",
            color="LightSkyBlue",
            label="Adjusted Incidence",
            markersize=8,
            markeredgewidth=1,
            markeredgecolor="white",
        )
        legend = ax2.legend(**plotting.LEGEND_OPTS)
        for handle in legend.legend_handles:
            handle.set_markersize(8)

        fig.tight_layout()
        return fig

    def parameter_df(self) -> pd.DataFrame:
        return pd.DataFrame(
            data=[(*k, *v) for k, v in self.adjustment_parameters.items()],
            columns="Species|Regression Method|A|b|sigma".split("|"),
        )

    def to_docx(
        self,
        report: Report | None = None,
        header_level: int = 1,
        show_title: bool = True,
    ):
        """Returns a word document report.

        Args:
            report (Report | None, optional): A optional report instance, otherwise create one.
            header_level (int, optional): The top-level header level, defaults to 1.
            show_title (bool, optional): Show the top level title, defaults True.
        """
        if report is None:
            report = Report.build_default()

        h1 = report.styles.get_header_style(header_level)
        h2 = report.styles.get_header_style(header_level)

        if show_title:
            report.document.add_paragraph("Rao-Scott Transformation", h1)

        report.document.add_paragraph("Summary", h2)
        write_setting_p(report, "Species: ", self.species.name.title())
        summary = self.summary_df()
        report.document.add_paragraph(df_to_table(report, summary))
        report.document.add_paragraph(
            add_mpl_figure(report.document, self.figure(figsize=(8, 4)), 6.5)
        )
        report.document.add_paragraph("Rao-Scott Transformation Parameters", h2)
        report.document.add_paragraph(
            df_to_table(report, self.parameter_df().query(f"Species == '{self.species}'"))
        )
        report.document.add_paragraph(
            "Fox JF, Hogan KA, Davis A. Dose-Response Modeling with Summary Data from Developmental Toxicity Studies. Risk Anal. 2017 May;37(5):905-917. PMID: 27567129. DOI: 10.1111/risa.12667."
        )
        report.document.add_paragraph("")

        # Check if any design effects were clipped to 1
        if (
            self.df["design_avg_below_1"].any()
            or self.df["design_ls_below_1"].any()
            or self.df["design_o_below_1"].any()
        ):
            footnote_text = (
                "Note: One or more dose groups had a design effect < 1. "
                "To prevent adjusted values greater than their original value, "
                "a lower limit of 1 was applied to the design effect."
            )
            report.document.add_paragraph(footnote_text)

        # Check summary table to see if any row in the Incidence column contains either an integer 0 or string "0"
        incidence_col = summary["Incidence"]
        has_zero = ((incidence_col == 0) | (incidence_col == "0")).any()
        if has_zero:
            footnote_text = (
                "Note: One or more dose groups have a fraction affected equal to zero; for these dose groups, no scaling was performed and the design effect was set equal to 1."
            )
            report.document.add_paragraph(footnote_text)

        return report.document

    def summary_df(self) -> pd.DataFrame:
        mapping = {
            "dose": "Dose",
            "n": "N",
            "incidence": "Incidence",
            "fraction_affected": "Fraction Affected",
            "design_ls": "Design Effect (LS)",
            "design_o": "Design Effect (OR)",
            "design_avg": "Design Effect (Average)",
            "n_adjusted": "N (Rao-Scott Transformed)",
            "incidence_adjusted": "Incidence (Rao-Scott Transformed)",
        }
        df = self.df[list(mapping.keys())].rename(columns=mapping)
        return df.fillna("--")

    def to_excel(self) -> BytesIO:
        """Returns an Excel report with worksheets summarizing the adjustment.

        Returns:
            BytesIO: An Excel worksheets.
        """
        f = BytesIO()
        with pd.ExcelWriter(f) as writer:
            for name, df in [("data", self.summary_df()), ("parameters", self.parameter_df())]:
                df.to_excel(writer, sheet_name=name, index=False)
        return f
