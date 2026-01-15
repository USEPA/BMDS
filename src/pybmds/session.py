from __future__ import annotations

import logging
from copy import copy, deepcopy
from itertools import cycle
from typing import Any, ClassVar

import numpy as np
import numpy.typing as npt
import pandas as pd

from . import __version__, bmdscore, constants, plotting
from .constants import (
    BMDS_BLANK_VALUE,
    MAXIMUM_POLYNOMIAL_ORDER,
    Dtype,
    Models,
    PriorClass,
    DistType,
)
from .datasets.base import DatasetSchemaBase, DatasetType
from .models import continuous as c3
from .models import dichotomous as d3
from .models import ma
from .models import cma
from .models import nested_dichotomous as nd3
from .models.base import BmdModel, BmdModelAveraging, BmdModelAveragingSchema, BmdModelSchema
from .recommender import Recommender, RecommenderSettings
from .reporting.styling import (
    Report,
    add_mpl_figure,
    df_to_table,
    plot_dr,
    write_base_frequentist_table,
    write_bayesian_table,
    write_citation,
    write_dataset_metadata,
    write_dataset_table,
    write_inputs_table,
    write_model,
    write_models,
)
from .selected import SelectedModel
from .types import session as schema

logger = logging.getLogger(__name__)


class Session:
    """A Session is bmd modeling session for a single dataset.

    The session contains the dataset, model configuration and results, and model recommendations
    and potentially model averaging results too. Sessions are a primary data type that
    should be able to be serialized and deserialized.
    """

    model_options: ClassVar = {
        Dtype.DICHOTOMOUS: {
            Models.Logistic: d3.Logistic,
            Models.LogLogistic: d3.LogLogistic,
            Models.Probit: d3.Probit,
            Models.LogProbit: d3.LogProbit,
            Models.QuantalLinear: d3.QuantalLinear,
            Models.Multistage: d3.Multistage,
            Models.Gamma: d3.Gamma,
            Models.Weibull: d3.Weibull,
            Models.DichotomousHill: d3.DichotomousHill,
        },
        Dtype.CONTINUOUS: {
            Models.Linear: c3.Linear,
            Models.Polynomial: c3.Polynomial,
            Models.Power: c3.Power,
            Models.Hill: c3.Hill,
            Models.ExponentialM3: c3.ExponentialM3,
            Models.ExponentialM5: c3.ExponentialM5,
        },
        Dtype.CONTINUOUS_INDIVIDUAL: {
            Models.Linear: c3.Linear,
            Models.Polynomial: c3.Polynomial,
            Models.Power: c3.Power,
            Models.Hill: c3.Hill,
            Models.ExponentialM3: c3.ExponentialM3,
            Models.ExponentialM5: c3.ExponentialM5,
        },
        Dtype.NESTED_DICHOTOMOUS: {
            Models.NestedLogistic: nd3.NestedLogistic,
            Models.NCTR: nd3.Nctr,
        },
    }

    def _ensure_continuous_ma_models(
        self,
        settings: dict | None = None,
        include_efsa: bool = False,
    ):
        """
        Ensure that all continuous models required for CMA exist in session.models.

        - BMDS models are always included:
          Power, Hill (CV + NCV)
          Exp3, Exp5 (CV + NCV + Lognormal)

        - EFSA models are included ONLY if include_efsa=True
        """
        settings = deepcopy(settings) if settings else {}

        suffix = {
            DistType.normal: "CV",
            DistType.normal_ncv: "NCV",
            DistType.log_normal: "Lognormal",
        }

        def _add(Model, disttypes):
            for dt in disttypes:
                # skip if already present
                for m in self.models:
                    if isinstance(m, Model) and m.settings.disttype == dt:
                        break
                else:
                    model_settings = deepcopy(settings)
                    model_settings["disttype"] = dt
                    model_settings.setdefault("name", f"{Model.__name__} ({suffix[dt]})")
                    self.models.append(Model(self.dataset, settings=model_settings))

        # --- BMDS models (always in CMA) ---
        _add(c3.Power, [DistType.normal, DistType.normal_ncv])
        _add(c3.Hill, [DistType.normal, DistType.normal_ncv])
        _add(
            c3.ExponentialM3,
            [DistType.normal, DistType.normal_ncv, DistType.log_normal],
        )
        _add(
            c3.ExponentialM5,
            [DistType.normal, DistType.normal_ncv, DistType.log_normal],
        )

        # --- EFSA models (opt-in only) ---
        if include_efsa:
            for Model in (
                c3.MultiplicativeHill,
                c3.InverseExponential,
                c3.Lognormal,
                c3.Gamma,
                c3.LMS,
            ):
                _add(
                    Model,
                    [DistType.normal, DistType.normal_ncv, DistType.log_normal],
                )

    def __init__(
        self,
        dataset: DatasetType,
        recommendation_settings: RecommenderSettings | None = None,
        id: int | None = None,
        name: str = "",
        description: str = "",
    ):
        self.id = id
        self.name = name
        self.description = description
        self.dataset = dataset
        self.models: list[BmdModel] = []
        self.ma_weights: npt.NDArray | None = None
        self.model_average: BmdModelAveraging | None = None
        self.recommendation_settings: RecommenderSettings | None = recommendation_settings
        self.recommender: Recommender | None = None
        self.selected: SelectedModel = SelectedModel(self)

    def add_default_bayesian_models(self, settings: dict | None = None, model_average: bool = True):
        settings = deepcopy(settings) if settings else {}
        settings["priors"] = PriorClass.bayesian

        if self.dataset.dtype in (Dtype.CONTINUOUS, Dtype.CONTINUOUS_INDIVIDUAL):
            # BMDS CMA defaults (with variance variants)
            self._add_bmds_continuous_ma_variants(settings=settings)

        else:
            for name in self.model_options[self.dataset.dtype].keys():
                model_settings = deepcopy(settings)
                if name in Models.VARIABLE_POLYNOMIAL():
                    model_settings.update(degree=2)
                self.add_model(name, settings=model_settings)

        if model_average and self.dataset.dtype in (
            constants.Dtype.DICHOTOMOUS,
            constants.Dtype.CONTINUOUS,
            constants.Dtype.CONTINUOUS_INDIVIDUAL,
        ):
            self.add_model_averaging()

    def add_default_models(self, settings: dict | None = None):
        for name in self.model_options[self.dataset.dtype].keys():
            model_settings = deepcopy(settings) if settings is not None else None
            if name in Models.VARIABLE_POLYNOMIAL():
                min_poly_order = 2
                max_poly_order = min(self.dataset.num_dose_groups - 1, MAXIMUM_POLYNOMIAL_ORDER + 1)
                for i in range(min_poly_order, max_poly_order):
                    poly_model_settings = (
                        deepcopy(model_settings) if model_settings is not None else {}
                    )
                    poly_model_settings["degree"] = i
                    self.add_model(name, settings=poly_model_settings)
            else:
                self.add_model(name, settings=model_settings)

    def add_model(self, name, settings=None):
        Model = self.model_options[self.dataset.dtype][name]
        instance = Model(dataset=self.dataset, settings=settings)
        self.models.append(instance)

    def set_ma_weights(self, weights: npt.ArrayLike | None = None):
        if weights is None:
            weights = np.full(len(self.models), 1 / len(self.models), dtype=np.float64)
        if len(self.models) != len(weights):
            raise ValueError(f"# model weights ({weights}) != num models {len(self.models)}")
        weights = np.array(weights)
        self.ma_weights = weights / weights.sum()

    def add_model_averaging(self, weights: list[float] | None = None, include_efsa: bool = False):
        """
        Must be added after other models are added since a shallow copy is taken, and the
        execution of model averaging assumes all other models were executed.
        """
        if weights or self.ma_weights is None:
            self.set_ma_weights(weights)

        if self.dataset.dtype is constants.Dtype.DICHOTOMOUS:
            instance = ma.BmdModelAveragingDichotomous(session=self, models=copy(self.models))

        elif self.dataset.dtype in (
            constants.Dtype.CONTINUOUS,
            constants.Dtype.CONTINUOUS_INDIVIDUAL,
        ):
            allowed_bmds = (
                c3.Power,
                c3.Hill,
                c3.ExponentialM3,
                c3.ExponentialM5,
            )

            allowed_efsa = (
                c3.MultiplicativeHill,
                c3.InverseExponential,
                c3.Lognormal,
                c3.Gamma,
                c3.LMS,
            )

            ma_models = [
                m
                for m in self.models
                if isinstance(m, allowed_bmds) or (include_efsa and isinstance(m, allowed_efsa))
            ]

            instance = cma.BmdModelAveragingContinuous(session=self, models=ma_models)

        else:
            raise ValueError(f"Model averaging not supported for dtype: {self.dataset.dtype}")

        self.model_average = instance

    def execute(self):
        # execute individual models
        for model in self.models:
            model.execute_job()

        # execute model average
        if self.model_average is not None:
            self.model_average.execute_job()

    @property
    def recommendation_enabled(self) -> bool:
        if self.recommender is None:
            self.recommender = Recommender(settings=self.recommendation_settings)
        return self.recommender.settings.enabled

    def recommend(self):
        if not self.recommendation_enabled or self.recommender is None:
            raise ValueError("Recommendation not enabled.")
        self.recommender.recommend(self.dataset, self.models)

    def select(self, model: BmdModel | None, notes: str = ""):
        self.selected.select(model, notes)

    @property
    def has_recommended_model(self) -> bool:
        return (
            self.recommendation_enabled
            and self.recommender.results.recommended_model_index is not None
        )

    @property
    def recommended_model(self) -> BmdModel | None:
        if self.has_recommended_model:
            return self.models[self.recommender.results.recommended_model_index]

    def accept_recommendation(self):
        """Select the recommended model, if one exists."""
        if self.has_recommended_model:
            index = self.recommender.results.recommended_model_index
            self.select(self.models[index], "Selected as best-fitting model")
        else:
            self.select(None, "No model was selected as a best-fitting model")

    def execute_and_recommend(self):
        self.execute()
        self.recommend()

    def is_bayesian(self) -> bool:
        """Determine if models are using a bayesian or frequentist approach."""
        if self.dataset.dtype == constants.Dtype.NESTED_DICHOTOMOUS:
            return False
        return self.models[0].settings.priors.is_bayesian

    def dll_version(self) -> str:
        return bmdscore.version()

    # serializing
    # -----------
    def serialize(self) -> SessionSchema:
        schema = SessionSchema(
            id=self.id,
            name=self.name,
            description=self.description,
            version=dict(
                python=__version__,
                dll=self.dll_version(),
            ),
            dataset=self.dataset.serialize(),
            models=[model.serialize() for model in self.models],
            selected=self.selected.serialize(),
        )
        if self.model_average is not None:
            schema.bmds_model_average = self.model_average.serialize(self)

        if self.recommender is not None:
            schema.recommender = self.recommender.serialize()

        return schema

    @classmethod
    def from_serialized(cls, data: dict) -> Session:
        try:
            dtype = data["dataset"]["dtype"]
        except KeyError as err:
            raise ValueError("Invalid JSON format") from err

        dataset = DatasetSchemaBase.get_subclass(dtype).model_validate(data["dataset"])
        model_base_class = BmdModelSchema.get_subclass(dtype)
        data["dataset"] = dataset
        data["models"] = [model_base_class.model_validate(model_) for model_ in data["models"]]
        ma = data.get("model_average")
        if ma:
            data["model_average"] = BmdModelAveragingSchema.get_subclass(dtype).model_validate(ma)
        return SessionSchema.model_validate(data).deserialize()

    # reporting
    # ---------
    def to_dict(self):
        return self.serialize().model_dump(by_alias=True)

    def session_title(self) -> str:
        if self.id and self.name:
            return f"${self.id}: {self.name}"
        elif self.name:
            return self.name
        elif self.id:
            return f"Session #{self.id}"
        elif self.dataset.metadata.name:
            return f"Session for {self.dataset.metadata.name}"
        else:
            return "Modeling Summary"

    def to_df(self, extras: dict | None = None, clean: bool = True) -> pd.DataFrame:
        """Export an executed session to a pandas dataframe.

        Args:
            extras (dict, optional): Extra items to add to row.
            clean (bool, default True): Remove empty columns.

        Returns:
            pd.DataFrame: A pandas dataframe
        """

        dataset_dict = {}
        self.dataset.update_record(dataset_dict)
        extras = extras or {}

        # add a row for each model
        models = []
        for bmds_model_index, model in enumerate(self.models):
            d: dict[str, Any] = {
                **extras,
                **dataset_dict,
                **dict(bmds_model_index=bmds_model_index, model_name=model.name()),
            }
            model.settings.update_record(d)
            model.results.update_record(d)

            if self.recommendation_enabled and self.recommender.results is not None:
                self.recommender.results.update_record(d, bmds_model_index)
                self.selected.update_record(d, bmds_model_index)

            if self.model_average:
                self.model_average.results.update_record_weights(d, bmds_model_index)

            models.append(d)

        # add model average row
        if self.model_average:
            d = dict(
                **extras,
                **dataset_dict,
                bmds_model_index=100,
                model_name="Model average",
            )
            self.model_average.settings.update_record(d)
            self.model_average.results.update_record(d)
            models.append(d)

        df = pd.DataFrame(data=models)
        if "slope_factor" in df.columns and np.allclose(df.slope_factor, BMDS_BLANK_VALUE):
            df.drop(columns=["slope_factor"], inplace=True)
        if clean:
            df = df.dropna(axis=1, how="all").fillna("")
        return df

    def to_docx(
        self,
        report: Report | None = None,
        header_level: int = 1,
        citation: bool = True,
        dataset_format_long: bool = True,
        all_models: bool = False,
        bmd_cdf_table: bool = False,
        session_inputs_table: bool = False,
    ):
        """Return a Document object with the session executed

        Args:
            report (Report, optional): A Report dataclass, or None to use default.
            header_level (int, optional): Starting header level. Defaults to 1.
            citation (bool, default True): Include citation
            dataset_format_long (bool, default True): long or wide dataset table format
            all_models (bool, default False):  Show all models, not just selected
            bmd_cdf_table (bool, default False): Export BMD CDF table
            session_inputs_table (bool, default False): Write an inputs table for a session,
                assuming a single model's input settings are representative of all models in a
                session, which may not always be true

        Returns:
            A python docx.Document object with content added.
        """
        if report is None:
            report = Report.build_default()

        h1 = report.styles.get_header_style(header_level)
        h2 = report.styles.get_header_style(header_level + 1)

        report.document.add_paragraph(self.session_title(), h1)
        if self.description:
            report.document.add_paragraph(self.description)

        report.document.add_paragraph("Dataset", h2)
        write_dataset_metadata(report, self.dataset)
        if len(self.models) > 0:
            first_model = self.models[0]
        write_dataset_table(report, self.dataset, long=dataset_format_long, model=first_model)

        if session_inputs_table:
            report.document.add_paragraph("Settings", h2)
            write_inputs_table(report, self)

        if self.is_bayesian():
            report.document.add_paragraph("Bayesian Summary", h2)
            write_bayesian_table(report, self)
            plot_dr(report, self)
            if self.model_average and bmd_cdf_table:
                report.document.add_paragraph("CDF:", report.styles.tbl_body)
                fig = self.model_average.cdf_plot(xlabel=self.dataset.get_xlabel())
                report.document.add_paragraph(add_mpl_figure(report.document, fig, 6))
                df_to_table(report, self.model_average.cdf())
            if all_models:
                report.document.add_paragraph("Individual Model Results", h2)
                write_models(report, self, bmd_cdf_table, header_level + 2)

        else:
            report.document.add_paragraph("Maximum Likelihood Approach", h2)
            write_base_frequentist_table(report, self)
            plot_dr(report, self)
            if all_models:
                report.document.add_paragraph("Individual Model Results", h2)
                write_models(report, self, bmd_cdf_table, header_level + 2)
            else:
                if self.selected.model:
                    write_model(
                        report,
                        self.selected.model,
                        bmd_cdf_table,
                        header_level + 1,
                        header_text=f"Selected Model: {self.selected.model.name()}",
                    )
                else:
                    report.document.add_paragraph("Selected Model", h2)
                    report.document.add_paragraph("No model was selected as a best-fitting model.")

        if citation:
            write_citation(report, header_level)

        return report.document

    def plot(self, figsize: tuple[float, float] | None = None, colorize: bool = True):
        """
        After model execution, print the dataset, curve-fit, BMD, and BMDL.
        """
        dataset = self.dataset
        fig = dataset.plot(figsize=figsize)
        ax = fig.gca()
        if self.dataset.dtype == constants.Dtype.DICHOTOMOUS:
            ax.set_ylim(-0.05, 1.05)
        has_ma = self.model_average is not None
        model_class = "MLE Models"
        if self.is_bayesian():
            model_class = "Bayesian Model Average" if has_ma else "Bayesian Models"
        title = f"{dataset._get_dataset_name()}\n{model_class}\n{self.models[0].settings.bmr_text}"
        ax.set_title(title)
        if colorize:
            color_cycle = cycle(plotting.INDIVIDUAL_MODEL_COLORS)
            line_cycle = cycle(plotting.INDIVIDUAL_LINE_STYLES)
        else:
            color_cycle = cycle(["#ababab"])
            line_cycle = cycle(["solid"])
        for i, model in enumerate(self.models):
            if colorize:
                label = model.name()
            elif i == 0:
                label = "Individual Model"
            else:
                label = None
            ax.plot(
                model.results.plotting.dr_x,
                model.results.plotting.dr_y,
                label=label,
                c=next(color_cycle),
                linestyle=next(line_cycle),
                zorder=100,
                lw=2,
            )
        if has_ma:
            ma = self.model_average
            ax.plot(
                ma.results.dr_x,
                ma.results.dr_y,
                label="Model Average",
                c="#6470C0",
                lw=4,
                zorder=110,
            )
            plotting.add_bmr_lines(
                ax, ma.results.bmd, ma.results.bmd_y, ma.results.bmdl, ma.results.bmdu
            )

        # reorder handles and labels
        handles, labels = ax.get_legend_handles_labels()
        if colorize:
            if "Fraction Affected ± 95% CI" in labels:
                idx = labels.index("Fraction Affected ± 95% CI")  # dichotomous
            elif "Observed Mean ± 95% CI" in labels:
                idx = labels.index("Observed Mean ± 95% CI")  # continuous summary
            elif "Observed" in labels:
                idx = labels.index("Observed")  # continuous individual
            elif "Fraction Affected" in labels:
                idx = labels.index("Fraction Affected")  # nested dichotomous
            else:
                raise ValueError("Cannot get label")
            labels.insert(0, labels.pop(idx))
            handles.insert(0, handles.pop(idx))

            if has_ma:
                idx = labels.index("Model Average")
                labels.append(labels.pop(idx))
                handles.append(handles.pop(idx))

                idx = labels.index("BMDL-BMD-BMDU")
                labels.append(labels.pop(idx))
                handles.append(handles.pop(idx))

            ax.legend(handles, labels, **plotting.LEGEND_OPTS)
        else:
            order = [2, 0, 1, 3] if has_ma else [1, 0]
            ax.legend(
                [handles[idx] for idx in order],
                [labels[idx] for idx in order],
                **plotting.LEGEND_OPTS,
            )
        if has_ma:
            plotting.improve_bmd_diamond_legend(ax)
        fig.tight_layout()
        return fig


class SessionSchema(schema.SessionSchemaBase):
    def deserialize(self) -> Session:
        session = Session(
            dataset=self.dataset.deserialize(),
            id=self.id,
            name=self.name,
            description=self.description,
        )
        session.models = [model.deserialize(session.dataset) for model in self.models]
        session.selected = self.selected.deserialize(session)
        if self.bmds_model_average is not None:
            session.model_average = self.bmds_model_average.deserialize(session)
        if self.recommender is not None:
            session.recommendation_settings = self.recommender.settings
            session.recommender = self.recommender.deserialize()
        return session
