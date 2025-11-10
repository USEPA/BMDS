from typing import Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel

#from .. import bmdscore
import bmdscore
from ..models.continuous import BmdModelContinuous
from .common import inspect_cpp_obj
from .continuous import NumpyFloatArray


class ModelAverageResult(BaseModel):
    pass


class ContinuousModelAverage:
    def __init__(self, dataset, models: list[BmdModelContinuous], model_weights: npt.NDArray):
        first = models[0].structs.analysis
        analysis = bmdscore.python_loud_continuous_analysis()
        analysis.BMD_type = first.BMD_type
        analysis.BMR = first.BMR
        analysis.alpha = first.alpha
        analysis.Y = dataset.incidences
        analysis.n_group = dataset.ns
        analysis.doses = dataset.doses
        analysis.n = dataset.num_dose_groups

        average = bmdscore.python_continuousMA_analysis()
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

        result = bmdscore.python_continuousMA_result()
        result.nmodels = len(models)
        result.dist_numE = 200
        result.models = [model.structs.result for model in models]
        result.bmdsRes = bmdsRes

        self.analysis = analysis
        self.average = average
        self.result = result
        self.bmdsRes = result.bmdsRes  # use this version; copied on assignment above

    def execute(self) -> "ContinuousModelAverageResult":
        bmdscore.pythonBMDSLoud(self.average, self.result)

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
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray

    @classmethod
    def from_cpp(cls, analysis: ContinuousModelAverage, model_results) -> Self:
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

if __name__ == "__main__":
    # Minimal dev test for the C++ dummy CMA hook.
    # This code runs only when you execute this file directly.

    import numpy as np

    # Import your dataset + session the same way you do in real code.
    # Adjust these imports to match your actual package layout.
    from .. import ContinuousDataset, Session  # or from ..session import Session, etc.

    from .continuous import NumpyFloatArray  # already imported above; safe if duplicated
    from .models.continuous import BmdModelContinuous  # adjust if needed, or reuse existing import

    # 1) Build a small continuous dataset (same as your example in the question)
    dataset = ContinuousDataset(
        doses=[0, 25, 50, 75, 100],
        ns=[20, 20, 20, 20, 20],
        means=[6, 8, 13, 25, 30],
        stdevs=[4, 4.3, 3.8, 4.4, 3.7],
    )

    # 2) Fit some continuous models the normal way
    session = Session(dataset=dataset)
    session.add_default_models()
    session.execute()

    # Keep only continuous models
    models = [m for m in session.models if isinstance(m, BmdModelContinuous)]

    # For now, just give them equal weights
    weights = np.full(len(models), 1 / len(models))

    # 3) Build the ContinuousModelAverage wrapper
    cma = ContinuousModelAverage(dataset, models, weights)

    # 4) This calls your C++ dummy:
    #    pythonBMDSLoud(self.average, self.result) -> pythonBMDSContLoud_dummy(...)
    cma.execute()

    # 5) Inspect what the C++ dummy wrote back
    print("=== Raw C++ dummy CMA values ===")
    print("BMD_MA:", cma.bmdsRes.BMD_MA)
    print("BMDL_MA:", cma.bmdsRes.BMDL_MA)
    print("BMDU_MA:", cma.bmdsRes.BMDU_MA)
    print("post_probs:", cma.result.post_probs)
    print("len(bmd_dist):", len(cma.result.bmd_dist))