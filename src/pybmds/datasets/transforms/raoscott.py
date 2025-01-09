"""
Rao Scott Adjustment, original citation:

Fox JF, Hogan KA, Davis A. Dose-Response Modeling with Summary Data from Developmental Toxicity
Studies. Risk Anal. 2017 May;37(5):905-917.
PMID: 27567129. DOI: 10.1111/risa.12667.
"""

from enum import StrEnum
from io import BytesIO
from typing import ClassVar, NamedTuple

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from ... import plotting
from ...reporting.styling import Report, add_mpl_figure, write_setting_p
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
        self.df = self.calculate()

    def calculate(self) -> pd.DataFrame:
        df = pd.DataFrame(
            {
                "dose": self.dataset.doses,
                "incidence": self.dataset.incidences,
                "n": self.dataset.ns,
            }
        )
        df["fraction_affected"] = df.incidence / df.n

        p = self.adjustment_parameters[(self.species, Regression.least_square)]
        df["design_ls"] = np.exp(p.a + (p.b * np.log(df.fraction_affected)) + (0.5 * p.sigma))

        p = self.adjustment_parameters[(self.species, Regression.orthogonal)]
        df["design_o"] = np.exp(p.a + (p.b * np.log(df.fraction_affected)) + (0.5 * p.sigma))

        df["design_avg"] = df[["design_ls", "design_o"]].mean(axis=1)
        df["incidence_adjusted"] = df.incidence / df.design_avg
        df["n_adjusted"] = df.n / df.design_avg
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
        h2 = report.styles.get_header_style(header_level + 1)

        if show_title:
            report.document.add_paragraph("Rao Scott Adjustment", h1)

        report.document.add_paragraph("Summary", h2)
        write_setting_p(report, "Species: ", self.species.name.title())

        report.document.add_paragraph("TODO - table")
        report.document.add_paragraph("TODO - figure")
        report.document.add_paragraph(
            add_mpl_figure(report.document, self.figure(figsize=(8, 4)), 6.5)
        )

        return report.document

    def to_excel(self) -> BytesIO:
        """Returns an Excel report with worksheets summarizing the adjustment.

        Returns:
            BytesIO: An Excel worksheets.
        """
        f = BytesIO()
        with pd.ExcelWriter(f) as writer:
            for name, df in [("data", self.df)]:
                df.to_excel(writer, sheet_name=name, index=False)
        return f
