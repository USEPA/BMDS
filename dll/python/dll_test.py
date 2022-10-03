#!/usr/bin/env python

import ctypes
import json
from enum import IntEnum
from pathlib import Path
from typing import Any, Dict, List, NamedTuple


def _list_to_c(list: List[Any], ctype):
    return (ctype * len(list))(*list)


class DichModel(IntEnum):
    d_hill = 1
    d_gamma = 2
    d_logistic = 3
    d_loglogistic = 4
    d_logprobit = 5
    d_multistage = 6
    d_probit = 7
    d_qlinear = 8
    d_weibull = 9


class DichotomousAnalysis(NamedTuple):
    model: int
    n: int
    Y: List[float]
    doses: List[float]
    n_group: List[float]
    prior: List[float]
    BMD_type: int
    BMR: float
    alpha: float
    degree: int
    samples: int
    burnin: int
    parms: int
    prior_cols: int

    class Struct(ctypes.Structure):
        _fields_ = [
            ("model", ctypes.c_int),  # Model Type as listed in DichModel
            ("n", ctypes.c_int),  # total number of observations obs/n
            ("Y", ctypes.POINTER(ctypes.c_double)),  # observed +
            ("doses", ctypes.POINTER(ctypes.c_double)),
            ("n_group", ctypes.POINTER(ctypes.c_double)),  # size of the group
            ("prior", ctypes.POINTER(ctypes.c_double)),  # a column order matrix parms X prior_cols
            ("BMD_type", ctypes.c_int),  # 1 = extra ; added otherwise
            ("BMR", ctypes.c_double),
            ("alpha", ctypes.c_double),  # alpha of the analysis
            ("degree", ctypes.c_int),  # degree of polynomial used only multistage
            ("samples", ctypes.c_int),  # number of MCMC samples
            ("burnin", ctypes.c_int),  # size of burnin
            ("parms", ctypes.c_int),  # number of parameters in the model
            ("prior_cols", ctypes.c_int),  # columns in the prior
        ]

        def dict(self) -> Dict:
            return dict(
                model=self.model,
                n=self.n,
                Y=self.Y[: self.n],
                doses=self.doses[: self.n],
                n_group=self.n_group[: self.n],
                prior=self.prior[: self.parms * self.prior_cols],
                BMD_type=self.BMD_type,
                BMR=self.BMR,
                alpha=self.alpha,
                degree=self.degree,
                samples=self.samples,
                burnin=self.burnin,
                parms=self.parms,
                prior_cols=self.prior_cols,
            )

    def to_c(self):
        return self.Struct(
            model=ctypes.c_int(self.model),
            n=ctypes.c_int(self.n),
            Y=_list_to_c(self.Y, ctypes.c_double),
            doses=_list_to_c(self.doses, ctypes.c_double),
            n_group=_list_to_c(self.n_group, ctypes.c_double),
            prior=_list_to_c(self.prior, ctypes.c_double),
            BMD_type=ctypes.c_int(self.BMD_type),
            BMR=ctypes.c_double(self.BMR),
            alpha=ctypes.c_double(self.alpha),
            degree=ctypes.c_int(self.degree),
            samples=ctypes.c_int(self.samples),
            burnin=ctypes.c_int(self.burnin),
            parms=ctypes.c_int(self.parms),
            prior_cols=ctypes.c_int(self.prior_cols),
        )


class DichotomousModelResult(NamedTuple):
    """
    Purpose: Data structure that is populated with all of the necessary
    information for a single model fit.
    """

    model: int
    nparms: int
    dist_numE: int

    class Struct(ctypes.Structure):

        _fields_ = [
            ("model", ctypes.c_int),  # dichotomous model specification
            ("nparms", ctypes.c_int),  # number of parameters in the model
            ("parms", ctypes.POINTER(ctypes.c_double)),  # parameter estimate
            ("cov", ctypes.POINTER(ctypes.c_double)),  # covariance estimate
            ("max", ctypes.c_double),  # value of the likelihood/posterior at the maximum
            ("dist_numE", ctypes.c_int),  # number of entries in rows for the bmd_dist
            ("model_df", ctypes.c_double),  # Used model degrees of freedom
            ("total_df", ctypes.c_double),  # Total degrees of freedom
            ("bmd_dist", ctypes.POINTER(ctypes.c_double)),  # bmd distribution (dist_numE x 2)
            ("bmd", ctypes.c_double),  # the central estimate of the BMD
        ]

        def dict(self) -> Dict:
            return dict(
                model=self.model,
                nparms=self.nparms,
                parms=self.parms[: self.nparms],
                cov=self.cov[: self.nparms ** 2],
                max=self.max,
                dist_numE=self.dist_numE,
                model_df=self.model_df,
                total_df=self.total_df,
                bmd_dist=self.bmd_dist[: self.dist_numE * 2],
            )

    def to_c(self):
        return self.Struct(
            model=ctypes.c_int(self.model),
            nparms=ctypes.c_int(self.nparms),
            parms=_list_to_c([0] * self.nparms, ctypes.c_double),
            cov=_list_to_c([0] * (self.nparms ** 2), ctypes.c_double),
            dist_numE=ctypes.c_int(self.dist_numE),
            bmd_dist=_list_to_c([0] * (self.dist_numE * 2), ctypes.c_double),
        )


def main():
    path = Path("/usr/local/lib/libDRBMD.so")
    assert path.exists()
    dll = ctypes.cdll.LoadLibrary(str(path))

    doses = [0, 50, 100, 150, 200]
    Y = [0, 5, 30, 65, 90]
    n_group = [100, 100, 100, 100, 100]
    prior = [1.0, 2.0, 0.0, 0.1, 2.0, 1.0, -20.0, 1e-12, 20.0, 100.0]
    prior_cols = 5
    parms = int(len(prior) / prior_cols)
    da = DichotomousAnalysis(
        model=DichModel.d_logistic.value,
        n=len(n_group),
        Y=Y,
        doses=doses,
        n_group=n_group,
        prior=prior,
        BMD_type=1,
        BMR=0.1,
        alpha=0.05,
        degree=parms - 1,
        samples=100,
        burnin=20,
        parms=parms,
        prior_cols=prior_cols,
    )
    da_res = DichotomousModelResult(
        model=DichModel.d_logistic.value, nparms=parms, dist_numE=200
    )

    da_struct = da.to_c()
    da_res_struct = da_res.to_c()
    dll.estimate_sm_laplace_dicho(
        ctypes.pointer(da_struct), ctypes.pointer(da_res_struct), True
    )

    print(
        json.dumps(
            dict(inputs=da_struct.dict(), outputs=da_res_struct.dict()), indent=2
        )
    )


if __name__ == "__main__":
    main()
