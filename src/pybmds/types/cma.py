from typing import Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel

from .. import bmdscore
from ..models.continuous import BmdModelContinuous
from .common import inspect_cpp_obj
from .continuous import NumpyFloatArray


class ModelAverageResult(BaseModel):
    pass


class ContinuousModelAverage:
    def __init__(
        self,
        dataset,
        models: list[BmdModelContinuous],
        model_weights: npt.NDArray,
        weight_option: int | None = None,
    ):
        first = models[0].structs.analysis
        analysis = bmdscore.python_continuous_analysis()
        analysis.BMD_type = first.BMD_type
        analysis.BMR = first.BMR
        analysis.alpha = first.alpha
        analysis.suff_stat = first.suff_stat
        analysis.disttype = first.disttype
        analysis.samples = first.samples
        analysis.burnin = first.burnin
        analysis.Y = dataset.means
        analysis.sd = dataset.stdevs
        analysis.n_group = dataset.ns
        analysis.doses = dataset.doses
        analysis.n = dataset.num_dose_groups

        average = bmdscore.python_continuousMA_analysis()
        average.nmodels = len(models)
        average.nparms = [model.structs.result.nparms for model in models]
        average.actual_parms = [model.structs.result.nparms for model in models]
        average.prior_cols = [model.structs.analysis.prior_cols for model in models]
        average.models = [int(m.structs.analysis.model) for m in models]
        average.disttype = [int(m.structs.analysis.disttype) for m in models]
        average.weightOption = 1 if weight_option is None else int(weight_option)

        if analysis.suff_stat is True:
            average.datatype = bmdscore.loud_datatype.l_summary
        elif analysis.suff_stat is False:
            average.datatype = bmdscore.loud_datatype.l_individual
        else:
            raise ValueError(f"Unsupported dataset dtype: {dataset.dtype}")

        average.priors = [model.structs.analysis.prior for model in models]

        average.modelPriors = [float(x) for x in model_weights]
        average.pyCA = analysis

        bmdsRes = bmdscore.BMDSMA_results()
        bmdsRes.BMD_MA = -9999.0
        bmdsRes.BMDL_MA = -9999.0
        bmdsRes.BMDU_MA = -9999.0
        bmdsRes.ebUpper = np.full(analysis.n, -9999)
        bmdsRes.ebLower = np.full(analysis.n, -9999)

        result = bmdscore.python_continuousMA_result()
        result.nmodels = len(models)
        result.dist_numE = 200

        result.post_probs = [0.0] * result.nmodels
        result.bmd_dist = [0.0] * int(analysis.samples)
        result.models = [bmdscore.python_continuous_model_result() for _ in range(result.nmodels)]

        for i, m in enumerate(models):
            nparms = int(m.structs.result.nparms)
            iters = int(analysis.samples)

            loud = result.models[i].loudRes

            loud.parms = [0.0] * (nparms * iters)

            loud.BMD = [0.0] * iters

        result.bmdsRes = bmdsRes

        self.analysis = analysis
        self.average = average
        self.result = result
        self.bmdsRes = result.bmdsRes  # use this version; copied on assignment above

    def execute(self) -> "ContinuousModelAverageResult":
        bmdscore.pythonBMDSLoud(self.average, self.result)
        return self

    def __str__(self) -> str:
        return "\n".join(
            [
                inspect_cpp_obj(self.analysis),
                inspect_cpp_obj(self.result),
            ]
        )


class ContinuousModelAverageResult(ModelAverageResult):
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
    model_bmd_dist: list[NumpyFloatArray]
    model_parm_dist: list[NumpyFloatArray]
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray

    @classmethod
    def from_cpp(cls, analysis: ContinuousModelAverage, model_results) -> Self:
        ma_bmd = np.asarray(analysis.result.bmd_dist, dtype=float)

        model_bmd: list[np.ndarray] = []
        model_parms: list[np.ndarray] = []

        for m in analysis.result.models:
            lr = m.loudRes  # fitResult

            b_raw = np.asarray(lr.BMD, dtype=float)
            n_draw = b_raw.size

            b = b_raw
            model_bmd.append(b)

            p = np.asarray(lr.parms, dtype=float)

            if p.ndim == 1 and n_draw > 0 and p.size % n_draw == 0:
                p = p.reshape(n_draw, p.size // n_draw)
            elif p.ndim == 2 and n_draw > 0 and p.shape[0] != n_draw and p.shape[1] == n_draw:
                p = p.T

            if p.ndim == 2:
                p = p[np.isfinite(p).all(axis=1)]
            else:
                p = p[np.isfinite(p)]

            model_parms.append(p)

        priors = np.asarray(analysis.average.modelPriors, dtype=float)
        posteriors = np.asarray(analysis.result.post_probs, dtype=float)
        if posteriors.size == 0:
            raise RuntimeError("C++ MA did not populate post_probs (empty).")

        values = np.asarray([r.plotting.dr_y for r in model_results], dtype=float)
        dr_x = np.asarray(model_results[0].plotting.dr_x, dtype=float)

        if posteriors.size != values.shape[0]:
            raise RuntimeError(f"post_probs ({posteriors.size}) != nmodels ({values.shape[0]}).")

        dr_y = posteriors @ values

        bmdsRes = analysis.result.bmdsRes
        bmds = np.asarray([bmdsRes.BMDL_MA, bmdsRes.BMD_MA, bmdsRes.BMDU_MA], dtype=float)
        bmds_ys = np.interp(bmds, dr_x, dr_y)

        return cls(
            bmdl=float(bmds[0]),
            bmd=float(bmds[1]),
            bmdu=float(bmds[2]),
            bmdl_y=float(bmds_ys[0]),
            bmd_y=float(bmds_ys[1]),
            bmdu_y=float(bmds_ys[2]),
            bmd_dist=ma_bmd,
            model_bmd_dist=model_bmd,
            model_parm_dist=model_parms,
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
