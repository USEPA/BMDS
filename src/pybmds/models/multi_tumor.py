from itertools import cycle
from typing import Self

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from pydantic import BaseModel

from .. import __version__, bmdscore, plotting, reporting
from ..constants import NUM_PRIOR_COLS
from ..datasets.dichotomous import DichotomousDataset, DichotomousDatasetSchema
from ..reporting.footnotes import TableFootnote
from ..reporting.styling import (
    Report,
    add_mpl_figure,
    set_column_width,
    write_bmd_cdf_table,
    write_cell,
    write_citation,
    write_pvalue_header,
)
from ..types.dichotomous import DichotomousModelSettings
from ..types.multi_tumor import MultitumorAnalysis, MultitumorResult, MultitumorSettings
from ..types.priors import multistage_cancer_prior
from ..types.session import VersionSchema
from ..utils import unique_items
from .dichotomous import MultistageCancer


def write_docx_frequentist_table(report: Report, session):
    """Add frequentist table to document."""
    styles = report.styles
    hdr = report.styles.tbl_header
    body = report.styles.tbl_body

    avg_row = len(session.datasets) > 1

    footnotes = TableFootnote()
    tbl = report.document.add_table(
        len(session.models) + 1 + (1 if avg_row else 0), 9, style=styles.table
    )

    write_cell(tbl.cell(0, 0), "Dataset: Model", style=hdr)
    write_cell(tbl.cell(0, 1), "BMDL", style=hdr)
    write_cell(tbl.cell(0, 2), "BMD", style=hdr)
    write_cell(tbl.cell(0, 3), "BMDU", style=hdr)
    write_cell(tbl.cell(0, 4), "CSF", style=hdr)
    write_pvalue_header(tbl.cell(0, 5), style=hdr)
    write_cell(tbl.cell(0, 6), "AIC", style=hdr)
    write_cell(tbl.cell(0, 7), "Scaled Residual at Control", style=hdr)
    write_cell(tbl.cell(0, 8), "Scaled Residual near BMD", style=hdr)

    if avg_row:
        write_cell(tbl.cell(1, 0), "Average", body)
        write_cell(tbl.cell(1, 1), session.results.bmdl, body)
        write_cell(tbl.cell(1, 2), session.results.bmd, body)
        write_cell(tbl.cell(1, 3), session.results.bmdu, body)
        write_cell(tbl.cell(1, 4), session.results.slope_factor, body)
        write_cell(tbl.cell(1, 5), "-", body)
        write_cell(tbl.cell(1, 6), "-", body)
        write_cell(tbl.cell(1, 7), "-", body)
        write_cell(tbl.cell(1, 8), "-", body)

    for ds_idx, ds_models in enumerate(session.models):
        row = ds_idx + 1 + (1 if avg_row else 0)
        idx = session.results.selected_model_indexes[ds_idx]
        if idx is None:
            dataset = session.datasets[ds_idx]
            write_cell(tbl.cell(row, 0), f"{dataset._get_dataset_name()}:\nno model", body)
            for j in range(1, 9):
                write_cell(tbl.cell(row, j), "-", body)
        else:
            model = ds_models[idx]
            write_cell(
                tbl.cell(row, 0), f"{model.dataset._get_dataset_name()}:\n{model.name()}", body
            )
            write_cell(tbl.cell(row, 1), model.results.bmdl, body)
            write_cell(tbl.cell(row, 2), model.results.bmd, body)
            write_cell(tbl.cell(row, 3), model.results.bmdu, body)
            write_cell(tbl.cell(row, 4), model.results.slope_factor, body)
            write_cell(tbl.cell(row, 5), model.get_gof_pvalue(), body)
            write_cell(tbl.cell(row, 6), model.results.fit.aic, body)
            write_cell(tbl.cell(row, 7), model.results.gof.residual[0], body)
            write_cell(tbl.cell(row, 8), model.results.gof.roi, body)

    # set column width
    widths = np.array([2, 1, 1, 1, 1, 1, 1, 1.5, 1.5])
    widths = widths / (widths.sum() / styles.portrait_width)
    for width, col in zip(widths, tbl.columns, strict=True):
        set_column_width(col, width)

    # write footnote
    if len(footnotes) > 0:
        footnotes.add_footnote_text(report.document, styles.tbl_footnote)


def write_docx_inputs_table(report: Report, session):
    """Add an input summary table to the document."""
    if len(session.models) == 0:
        raise ValueError("No models available")

    styles = report.styles
    hdr = report.styles.tbl_header
    body = report.styles.tbl_body

    settings = [models[len(models) - 1].settings for models in session.models]
    rows = {
        "Setting": "Value",
        "BMR": unique_items(settings, "bmr_text"),
        "Confidence Level (one sided)": unique_items(settings, "confidence_level"),
        "Maximum Degree": unique_items(settings, "degree"),
    }
    tbl = report.document.add_table(len(rows), 2, style=styles.table)
    for idx, (key, value) in enumerate(rows.items()):
        write_cell(tbl.cell(idx, 0), key, style=hdr)
        write_cell(tbl.cell(idx, 1), value, style=hdr if idx == 0 else body)


def write_docx_model(
    report: Report,
    model,
    header_level: int = 1,
    bmd_cdf_table: bool = False,
    include_dataset_name: bool = False,
):
    styles = report.styles
    header_style = styles.get_header_style(header_level)
    name = (
        f"{model.dataset._get_dataset_name()}: {model.name()}"
        if include_dataset_name
        else model.name()
    )
    report.document.add_paragraph(name, header_style)
    if model.has_results:
        report.document.add_paragraph(add_mpl_figure(report.document, model.plot(), 6))
        if bmd_cdf_table:
            report.document.add_paragraph(add_mpl_figure(report.document, model.cdf_plot(), 6))
        report.document.add_paragraph(model.text(), styles.fixed_width)
        if bmd_cdf_table:
            report.document.add_paragraph("CDF:", styles.tbl_body)
            write_bmd_cdf_table(report, model.results.fit)


class Multitumor:
    def __init__(
        self,
        datasets: list[DichotomousDataset],
        degrees: list[int] | None = None,
        settings: DichotomousModelSettings | dict | None = None,
        id: int | None = None,
        name: str = "",
        description: str = "",
        results: MultitumorResult | None = None,
    ):
        if len(datasets) == 0:
            raise ValueError("Must provide at least one dataset")
        self.id = id
        self.name = name
        self.description = description
        self.datasets = datasets
        for i, dataset in enumerate(datasets, start=1):
            if dataset.metadata.id is None:
                dataset.metadata.id = i
        self.degrees: list[int] = degrees or [0] * len(datasets)
        self.settings: DichotomousModelSettings = self.get_base_settings(settings)
        self.results = results
        self.structs: tuple | None = None
        self.models: list[list[MultistageCancer]] = []

    def get_base_settings(
        self, settings: DichotomousModelSettings | dict | None
    ) -> DichotomousModelSettings:
        if settings is None:
            return DichotomousModelSettings()
        elif isinstance(settings, DichotomousModelSettings):
            return settings
        else:
            return DichotomousModelSettings.model_validate(settings)

    def _build_model_settings(self) -> list[list[DichotomousModelSettings]]:
        # Build individual model settings based from inputs
        settings = []
        for i, dataset in enumerate(self.datasets):
            ds_settings = []
            degree_i = self.degrees[i]
            degrees_i = (
                range(degree_i, degree_i + 1)
                if degree_i > 0
                else range(1, min(dataset.num_dose_groups, 9))  # max of 8 if degree is 0 (auto)
            )
            for degree in degrees_i:
                model_settings = self.settings.model_copy(
                    update=dict(degree=degree, priors=multistage_cancer_prior())
                )
                ds_settings.append(model_settings)
            settings.append(ds_settings)
        return settings

    def to_cpp(self) -> MultitumorAnalysis:
        all_settings = self._build_model_settings()
        dataset_models = []
        dataset_results = []
        ns = []
        for dataset, dataset_settings in zip(self.datasets, all_settings, strict=True):
            mc_models = []
            self.models.append(mc_models)
            models = []
            results = []
            ns.append(dataset.num_dose_groups)
            for settings in dataset_settings:
                model = MultistageCancer(dataset, settings=settings)
                inputs = model._build_inputs()
                structs = inputs.to_cpp()
                models.append(structs.analysis)
                results.append(structs.result)
                mc_models.append(model)
            dataset_models.append(models)
            dataset_results.append(results)

        analysis = bmdscore.python_multitumor_analysis()
        analysis.BMD_type = self.settings.bmr_type.value
        analysis.BMR = self.settings.bmr
        analysis.alpha = self.settings.alpha
        analysis.degree = self.degrees
        analysis.models = dataset_models
        analysis.n = ns
        analysis.ndatasets = len(self.datasets)
        analysis.nmodels = [len(models) for models in dataset_models]
        analysis.prior_cols = NUM_PRIOR_COLS

        result = bmdscore.python_multitumor_result()
        result.ndatasets = len(self.datasets)
        result.nmodels = [len(results) for results in dataset_results]
        result.models = dataset_results

        return MultitumorAnalysis(analysis, result)

    def execute(self):
        self.structs = self.to_cpp()
        self.structs.execute()
        for i, models in enumerate(self.structs.result.models):
            for j, model in enumerate(models):
                if model.bmdsRes.validResult:
                    bmr = self.structs.analysis.models[i][j].BMR
                    model.bmdsRes.setSlopeFactor(bmr)
        self.results = MultitumorResult.from_model(self)
        return self.results

    def text(self) -> str:
        return self.results.text(self.datasets, self.models)

    def to_dict(self):
        return self.serialize().model_dump(by_alias=True)

    def serialize(self) -> "MultitumorSchema":
        return MultitumorSchema(
            version=self._serialize_version(),
            datasets=[ds.serialize() for ds in self.datasets],
            id=self.id,
            name=self.name,
            description=self.description,
            settings=self._serialize_settings(),
            results=self.results,
        )

    @classmethod
    def from_serialized(cls, data: dict) -> Self:
        return MultitumorSchema.model_validate(data).deserialize()

    def _serialize_version(self) -> VersionSchema:
        return VersionSchema(
            python=__version__,
            dll=bmdscore.version(),
        )

    def _serialize_settings(self) -> MultitumorSettings:
        return MultitumorSettings(
            degrees=self.degrees,
            bmr=self.settings.bmr,
            bmr_type=self.settings.bmr_type,
            alpha=self.settings.alpha,
        )

    def to_df(self, extras: dict | None = None, clean: bool = True) -> pd.DataFrame:
        """Export an executed session to a pandas dataframe.

        Args:
            extras (dict, optional): Extra items to add to row.
            clean (bool, default True): Remove empty columns.

        Returns:
            pd.DataFrame: A pandas dataframe
        """
        if extras is None:
            extras = {}
        results = self.results
        data = []

        # model average
        ma = extras.copy()
        ma.update(
            dataset_index=np.nan,
            dataset_id=np.nan,
            dataset_name=np.nan,
            dataset_dose_name=np.nan,
            dataset_dose_units=np.nan,
            dataset_response_name=np.nan,
            dataset_response_units=np.nan,
            dataset_doses=np.nan,
            dataset_ns=np.nan,
            dataset_incidences=np.nan,
            model_index=np.nan,
            model_name="Model average",
            slope_factor=results.slope_factor,
            selected=np.nan,
            bmdl=results.bmdl,
            bmd=results.bmd,
            bmdu=results.bmdu,
            aic=np.nan,
            loglikelihood=np.nan,
            p_value=np.nan,
            overall_dof=np.nan,
            bic_equiv=np.nan,
            chi_squared=np.nan,
            residual_of_interest=np.nan,
            residual_at_lowest_dose=np.nan,
        )
        data.append(ma)

        # add models
        for dataset_i, models in enumerate(results.models):
            dataset = self.datasets[dataset_i]
            extras.update(dataset_index=dataset_i)
            dataset.update_record(extras)
            # individual model rows
            for model_i, model in enumerate(models):
                extras.update(
                    model_index=model_i,
                    model_name=model.name,
                    selected=results.selected_model_indexes[dataset_i] == model_i,
                )
                model.results.update_record(extras)
                data.append(extras.copy())

        df = pd.DataFrame(data=data)
        if clean:
            df = df.dropna(axis=1, how="all").fillna("")
        return df

    def params_df(self, extras: dict | None = None) -> pd.DataFrame:
        """Returns a pd.DataFrame of all parameters for all models executed.

        Args:
            extras (dict | None): extra columns to prepend
        """
        data = []
        extras = extras or {}
        for dataset_index, dataset_models in enumerate(self.results.models):
            dataset = self.datasets[dataset_index]
            for model_index, model in enumerate(dataset_models):
                data.extend(
                    model.results.parameters.rows(
                        extras={
                            **extras,
                            "dataset_id": dataset.metadata.id,
                            "dataset_name": dataset.metadata.name,
                            "model_index": model_index,
                            "model_name": model.name,
                        }
                    )
                )
        return pd.DataFrame(data)

    def datasets_df(self, extras: dict | None = None) -> pd.DataFrame:
        """Returns a pd.DataFrame of all datasets within a session.

        Args:
            extras (dict | None): extra columns to prepend
        """

        data = []
        for dataset in self.datasets:
            data.extend(dataset.rows(extras))
        return pd.DataFrame(data)

    def session_title(self) -> str:
        if self.id and self.name:
            return f"${self.id}: {self.name}"
        elif self.name:
            return self.name
        elif self.id:
            return f"Session #{self.id}"
        else:
            return "Modeling Summary"

    def dll_version(self) -> str:
        return bmdscore.version()

    def plot(self, figsize: tuple[float, float] | None = None) -> Figure:
        fig = plotting.create_empty_figure(figsize=figsize)
        ax = fig.axes[0]

        # add individual model fits
        selected_models = []
        color_cycle = cycle(plotting.INDIVIDUAL_MODEL_COLORS)
        for idx, dataset in enumerate(self.datasets):
            color = next(color_cycle)
            if idx == 0:
                ax.set_xlabel(dataset.get_xlabel())
                ax.set_ylabel(dataset.get_ylabel())
                ax.set_title("Model Average")
            selected_idx = self.results.selected_model_indexes[idx]
            ax.scatter(
                dataset.doses,
                dataset.plot_data().mean,
                c="white",
                edgecolors=color,
                s=70,
                linewidth=2,
                label=dataset._get_dataset_name() if selected_idx is None else None,
            )
            if selected_idx is not None:
                model = self.models[idx][selected_idx]
                ax.plot(
                    model.results.plotting.dr_x,
                    model.results.plotting.dr_y,
                    label=f"{dataset._get_dataset_name()}; {model.name()}",
                    c=color,
                )
                selected_models.append(model)

        # add slope factor line and bmd interval
        if self.results.bmdl and self.results.slope_factor:
            x0 = np.mean([model.results.plotting.dr_x[0] for model in selected_models])
            y0 = np.mean([model.results.plotting.dr_y[0] for model in selected_models])
            y1 = y0 + self.results.bmdl * self.results.slope_factor
            plotting.add_bmr_lines(
                ax,
                self.results.bmd,
                y1,
                self.results.bmdl,
                self.results.bmdu,
                c="k",
                ecolor=plotting.to_rgba("k", 0.7),
            )
            ax.plot(
                [x0, self.results.bmdl],
                [y0, y1],
                label="Slope Factor",
                **{**plotting.LINE_FORMAT, "linestyle": (0, (2, 1, 1, 1)), "c": "#ef7215"},
            )

        ax.legend(**plotting.LEGEND_OPTS)
        plotting.improve_bmd_diamond_legend(ax)

        return fig

    def to_docx(
        self,
        report: Report | None = None,
        header_level: int = 1,
        citation: bool = False,
        dataset_format_long: bool = True,
        all_models: bool = False,
        bmd_cdf_table: bool = False,
        **kw,
    ):
        """Return a Document object with the session executed

        Args:
            report (Report, optional): A Report dataclass, or None to use default.
            header_level (int, optional): Starting header level. Defaults to 1.
            citation (bool, default True): Include citation
            dataset_format_long (bool, default True): long or wide dataset table format
            all_models (bool, default False):  Show all models, not just selected
            bmd_cdf_table (bool, default False): Export BMD CDF table

        Returns:
            A python docx.Document object with content added
        """
        if report is None:
            report = Report.build_default()

        h1 = report.styles.get_header_style(header_level)
        h2 = report.styles.get_header_style(header_level + 1)
        h3 = report.styles.get_header_style(header_level + 2)
        report.document.add_paragraph(self.session_title(), h1)

        if self.description:
            report.document.add_paragraph(self.description)

        report.document.add_paragraph("Datasets", h2)
        for dataset in self.datasets:
            reporting.styling.write_dataset_metadata(report, dataset)
            reporting.styling.write_dataset_table(report, dataset, dataset_format_long)

        report.document.add_paragraph("Settings", h2)
        write_docx_inputs_table(report, self)

        report.document.add_paragraph("Maximum Likelihood Approach", h2)
        write_docx_frequentist_table(report, self)
        report.document.add_paragraph(add_mpl_figure(report.document, self.plot(), 6))
        report.document.add_paragraph(self.results.ms_combo_text(), report.styles.fixed_width)

        report.document.add_paragraph("Individual Model Results", h2)

        for dataset, selected_idx, models in zip(
            self.datasets, self.results.selected_model_indexes, self.models, strict=True
        ):
            if all_models:
                report.document.add_paragraph(dataset._get_dataset_name(), h3)
                for model in models:
                    write_docx_model(
                        report,
                        model,
                        header_level=header_level + 2,
                        bmd_cdf_table=bmd_cdf_table,
                        include_dataset_name=False,
                    )
            else:
                if selected_idx is None:
                    report.document.add_paragraph(dataset._get_dataset_name(), h3)
                    report.document.add_paragraph("No model was selected.")
                else:
                    model = models[selected_idx]
                    write_docx_model(
                        report,
                        model,
                        header_level=header_level + 2,
                        bmd_cdf_table=bmd_cdf_table,
                        include_dataset_name=True,
                    )

        if citation:
            write_citation(report, header_level=header_level)

        return report.document


class MultitumorSchema(BaseModel):
    version: VersionSchema
    id: int | None
    name: str = ""
    description: str = ""
    datasets: list[DichotomousDatasetSchema]
    settings: MultitumorSettings
    results: MultitumorResult | None

    def deserialize(self) -> Multitumor:
        datasets = [ds.deserialize() for ds in self.datasets]
        settings = dict(
            bmr=self.settings.bmr,
            bmr_type=self.settings.bmr_type,
            alpha=self.settings.alpha,
        )
        mt = Multitumor(
            datasets=datasets,
            degrees=self.settings.degrees,
            settings=settings,
            id=self.id,
            name=self.name,
            description=self.description,
            results=self.results,
        )
        # hydrate models
        for dataset, ds_models in zip(mt.datasets, mt.results.models, strict=True):
            mt.models.append([MultistageCancer.deserialize(dataset, model) for model in ds_models])
        return mt
