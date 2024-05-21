from typing import Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel

from .. import bmdscore
from ..models.dichotomous import BmdModelDichotomous
from .common import inspect_cpp_obj
from .continuous import NumpyFloatArray


class ModelAverageResult(BaseModel):
    pass


class DichotomousModelAverage:
    def __init__(self, dataset, models: list[BmdModelDichotomous], model_weights: npt.NDArray):
        first = models[0].structs.analysis
        analysis = bmdscore.python_dichotomous_analysis()
        analysis.BMD_type = first.BMD_type
        analysis.BMR = first.BMR
        analysis.alpha = first.alpha
        analysis.Y = dataset.incidences
        analysis.n_group = dataset.ns
        analysis.doses = dataset.doses
        analysis.n = dataset.num_dose_groups

        average = bmdscore.python_dichotomousMA_analysis()
        average.nmodels = len(models)
        average.nparms = [model.structs.result.nparms for model in models]
        average.actual_parms = [model.structs.result.nparms for model in models]
        average.prior_cols = [model.structs.analysis.prior_cols for model in models]
        average.models = [model.structs.analysis.model for model in models]
        average.priors = [model.structs.analysis.prior for model in models]
        average.modelPriors = model_weights
        average.pyDA = analysis

        bmdsRes = bmdscore.BMDSMA_results()
        bmdsRes.BMD = np.full(average.nmodels, -9999)
        bmdsRes.BMDL = np.full(average.nmodels, -9999)
        bmdsRes.BMDU = np.full(average.nmodels, -9999)
        bmdsRes.ebUpper = np.full(analysis.n, -9999)
        bmdsRes.ebLower = np.full(analysis.n, -9999)

        result = bmdscore.python_dichotomousMA_result()
        result.nmodels = len(models)
        result.dist_numE = 200
        result.models = [model.structs.result for model in models]
        result.bmdsRes = bmdsRes

        self.analysis = analysis
        self.average = average
        self.result = result
        self.bmdsRes = result.bmdsRes  # use this version; copied on assignment above

    def execute(self) -> "DichotomousModelAverageResult":
        bmdscore.pythonBMDSDichoMA(self.average, self.result)

    def __str__(self) -> str:
        lines = []
        inspect_cpp_obj(lines, self.analysis, depth=0)
        inspect_cpp_obj(lines, self.result, depth=0)
        return "\n".join(lines)


class DichotomousModelAverageResult(ModelAverageResult):
    """
    Model average fit
    """

    bmdl: float
    bmd: float
    bmdu: float
    bmdl_y: float
    bmd_y: float
    bmdu_y: float
    bmd_dist: NumpyFloatArray
    priors: NumpyFloatArray
    posteriors: NumpyFloatArray
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray

    @classmethod
    def from_cpp(cls, analysis: DichotomousModelAverage, model_results) -> Self:
        # only keep positive finite values
        arr = np.array(analysis.result.bmd_dist).reshape(2, analysis.result.dist_numE).T
        arr = arr[np.isfinite(arr[:, 0])]
        arr = arr[arr[:, 0] > 0]

        # calculate dr_y for model averaging
        priors = np.array(analysis.average.modelPriors)
        posteriors = np.array(analysis.result.post_probs)
        values = np.array([result.plotting.dr_y for result in model_results])
        dr_x = model_results[0].plotting.dr_x
        dr_y = values.T.dot(posteriors)
        bmds = [analysis.bmdsRes.BMDL_MA, analysis.bmdsRes.BMD_MA, analysis.bmdsRes.BMDU_MA]
        bmds_ys = np.interp(bmds, dr_x, dr_y)
        return cls(
            bmdl=bmds[0],
            bmd=bmds[1],
            bmdu=bmds[2],
            bmdl_y=bmds_ys[0],
            bmd_y=bmds_ys[1],
            bmdu_y=bmds_ys[2],
            bmd_dist=arr.T,
            priors=priors,
            posteriors=posteriors,
            dr_x=dr_x,
            dr_y=dr_y,
        )

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            bmdl=self.bmdl,
            bmd=self.bmd,
            bmdu=self.bmdu,
        )

    def update_record_weights(self, d: dict, index: int) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            model_prior=self.priors[index],
            model_posterior=self.posteriors[index],
        )
