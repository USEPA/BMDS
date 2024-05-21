"""
Poly K adjustment test, original citation:

Bailer AJ, Portier CJ. Effects of treatment-induced mortality and tumor-induced mortality on tests
for carcinogenicity in small samples. Biometrics. 1988 Jun;44(2):417-31.
PMID: 3390507. DOI: 10.2307/2531856
"""

from io import BytesIO
from itertools import cycle

import matplotlib.ticker as mtick
import pandas as pd
from matplotlib.figure import Figure

from ... import plotting
from ...reporting.styling import Report, add_mpl_figure, df_to_table, write_setting_p


class PolyKAdjustment:
    def __init__(
        self,
        doses: list[float],
        day: list[int],
        has_tumor: list[int],
        k: float = 3,
        max_day: float | None = None,
        dose_units: str = "",
        dose_name: str = "Dose",
        time_units: str = "days",
    ) -> None:
        """Adjust tumor data via poly K adjustment.

        Args:
            doses (list[float]): a list of floats
            day (list[int]): a list of days or other time measurements
            has_tumor (list[int]): integer, 0 or 1, where 1 is has_tumor is True
            k (float | None, optional): The adjustment term to apply. Defaults to 3.
            max_day (float | None, optional): The maximum data. If specific, the value is used,
                otherwise, it is calculated from the maximum reported day in the dataset.
            dose_units (str): Dose Units, for graphing
            dose_name (str): the term for the dose axis. Defaults to "Dose".
            time_units (str): the time units, defaults to "days".
        """
        self.k = k
        self.max_day = max_day
        self.dose_units = dose_units
        self.dose_name = dose_name
        self.time_units = time_units
        self.input_data = pd.DataFrame(dict(dose=doses, day=day, has_tumor=has_tumor))
        self.adjusted_data = self.calc_adjusted_n()
        self.summary = self.calc_summary_stats()

    def calc_adjusted_n(self) -> pd.DataFrame:
        """Adjust N for individual observations in a dataset."""
        columns = ["dose", "day", "has_tumor"]
        if self.input_data.columns.tolist() != columns:
            raise ValueError(f"Unexpected column names; expecting {columns}")
        if set(self.input_data.has_tumor.unique()) != {0, 1}:
            raise ValueError("Expected `has_tumor` values must be 0 and 1")
        df = self.input_data.copy()
        max_day = self.max_day or df.day.max()
        df.loc[:, "adj_n"] = (df.query("has_tumor==0").day / max_day) ** self.k
        df.loc[:, "adj_n"] = df.loc[:, "adj_n"].fillna(1).clip(upper=1)
        return df

    def calc_summary_stats(self) -> pd.DataFrame:
        """Calculate grouped summary statistics by group for Poly K adjusted data."""
        columns = ["dose", "day", "has_tumor", "adj_n"]
        if self.adjusted_data.columns.tolist() != columns:
            raise ValueError(f"Unexpected column names: {columns}")
        grouped = self.adjusted_data.groupby("dose")
        df2 = pd.DataFrame(
            data=[
                grouped.has_tumor.count().rename("n", inplace=True),
                grouped.adj_n.sum().rename("adj_n", inplace=True),
                grouped.has_tumor.sum().rename("incidence", inplace=True),
            ]
        ).T.reset_index()
        df2.loc[:, "proportion"] = df2.incidence / df2.n
        df2.loc[:, "adj_proportion"] = df2.incidence / df2.adj_n
        return df2

    def summary_figure(self, figsize: tuple[float, float] | None = None) -> Figure:
        fig = plotting.create_empty_figure(figsize=figsize)
        ax = fig.gca()

        ax.set_title("Adjusted Proportion vs Original Proportion")
        ax.set_xlabel(
            f"{self.dose_name} ({self.dose_units})" if self.dose_units else self.dose_name
        )
        ax.set_ylabel("Proportion (%)")
        ax.margins(plotting.PLOT_MARGINS)

        ax.plot(
            self.summary.dose,
            self.summary.proportion,
            "o-",
            color="blue",
            label="Original Proportion",
            markersize=8,
            markeredgewidth=1,
            markeredgecolor="white",
        )
        ax.plot(
            self.summary.dose,
            self.summary.adj_proportion,
            "^-",
            color="red",
            label="Adjusted Proportion",
            markersize=8,
            markeredgewidth=1,
            markeredgecolor="white",
        )

        legend = ax.legend(**plotting.LEGEND_OPTS)
        for handle in legend.legend_handles:
            handle.set_markersize(8)

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
        return fig

    def tumor_incidence_figure(self, figsize: tuple[float, float] | None = None) -> Figure:
        markers = "o^svDP"
        df = self.input_data.copy()
        df.loc[:, "cummulative_tumor"] = df.groupby("dose").has_tumor.cumsum()

        marker = cycle(markers)

        fig = plotting.create_empty_figure(figsize=figsize)
        ax = fig.gca()
        ax.set_xlabel(f"Study duration ({self.time_units})")
        ax.set_ylabel("Cumulative tumor incidence")
        ax.margins(plotting.PLOT_MARGINS)
        ax.set_title("Tumor incidence over study duration")

        dose_units = f" {self.dose_units}" if self.dose_units else ""

        for value, d in df.groupby("dose"):
            ax.plot(
                d.day,
                d.cummulative_tumor,
                f"{next(marker)}-",
                label=str(value) + dose_units,
                markersize=8,
                markeredgewidth=1,
                markeredgecolor="white",
            )

        legend = ax.legend(**plotting.LEGEND_OPTS)
        for handle in legend.legend_handles:
            handle.set_markersize(8)

        return fig

    def write_docx_adjustment_table(self, report: Report):
        """Add adjusted input data table to the document."""
        df_to_table(report, self.adjusted_data)

    def write_docx_summary_table(self, report: Report):
        """Add a 'result'' data table with adjusted figures to the document."""
        df_to_table(report, self.summary)

    def to_docx(
        self,
        report: Report | None = None,
        header_level: int = 1,
        show_title: bool = True,
    ):
        """Returns a word document report of the Poly K Adjustment Calculation.

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
            report.document.add_paragraph("Poly K Adjustment", h1)

        report.document.add_paragraph("Summary", h2)
        write_setting_p(report, "K Adjustment Factor: ", str(self.k))
        if self.max_day:
            write_setting_p(report, "Maximum Day: ", str(self.max_day))

        report.document.add_paragraph(self.write_docx_summary_table(report))
        report.document.add_paragraph(add_mpl_figure(report.document, self.summary_figure(), 6))

        report.document.add_paragraph("Individual Adjustments", h2)
        report.document.add_paragraph(
            add_mpl_figure(report.document, self.tumor_incidence_figure(), 6)
        )
        self.write_docx_adjustment_table(report)

        return report.document

    def to_excel(self) -> BytesIO:
        """Returns an Excel report with worksheets summarizing the adjustment.

        Returns:
            BytesIO: An Excel worksheets.
        """
        f = BytesIO()
        with pd.ExcelWriter(f) as writer:
            for name, df in [("adjusted", self.adjusted_data), ("summary", self.summary)]:
                df.to_excel(writer, sheet_name=name, index=False)
        return f
