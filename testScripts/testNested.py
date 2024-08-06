import sys

import bmdscore
import numpy as np


def getNLogisticPrior(ngrp, prior_cols, restricted):
    prG = np.array([0, 1])
    prB = np.array([-1 * sys.float_info.max, sys.float_info.max])
    prT1 = np.array([0, 1])
    prT2 = np.array([-1 * sys.float_info.max, sys.float_info.max])
    prUR = np.array([0, 18])
    prRR = np.array([1, 18])
    prP = np.array([0, sys.float_info.max])
    pr = []
    for i in range(prior_cols):
        pr.append(prG[i])
        pr.append(prB[i])
        pr.append(prT1[i])
        pr.append(prT2[i])
        if restricted:
            pr.append(prRR[i])
        else:
            pr.append(prUR[i])
        for j in range(ngrp):
            pr.append(prP[i])
    return pr


pyAnal = bmdscore.python_nested_analysis()
pyAnal.model = bmdscore.nested_model.nlogistic
pyAnal.restricted = True

# fmt: off
pyAnal.doses = np.array(
    [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        100, 100, 100, 100, 100, 100, 100, 100, 100
    ]
)
pyAnal.litterSize = np.array(
    [
        16, 9, 15, 14, 13, 9, 10, 14, 10, 11, 14, 9, 14, 9, 13, 12, 10, 10, 11, 14,
        11, 11, 14, 11, 10, 11, 10, 15, 7, 14, 11, 14, 12, 13, 12, 14, 11, 8, 10
    ]
)
pyAnal.incidence = np.array(
    [
        1, 1, 2, 3, 3, 0, 2, 2, 1, 2, 4, 5, 6, 2, 6, 3, 1, 2, 4,
        3, 4, 5, 5, 4, 5, 4, 5, 6, 2, 4, 6, 6, 8, 7, 8, 6, 6, 5, 4
    ]
)
pyAnal.lsc = np.array(
    [
        16, 9, 15, 14, 13, 9, 10, 14, 10, 11, 14, 9, 14, 9, 13, 12, 10, 10,
        11, 14, 11, 11, 14, 11, 10, 11, 10, 15, 7, 14, 11, 14, 12, 13, 12, 14, 11, 8, 10,
    ]
)
# fmt: on

pyAnal.LSC_type = 1
pyAnal.ILC_type = 1
pyAnal.BMD_type = 1
pyAnal.estBackground = True
pyAnal.BMR = 0.1
pyAnal.alpha = 0.05
pyAnal.numBootRuns = 3
pyAnal.iterations = 1000
pyAnal.seed = -9999
pyAnal.prior_cols = 2

ngrp = len(np.unique(pyAnal.doses))
Nobs = len(pyAnal.doses)

pyAnal.parms = 5 + ngrp
pyAnal.prior = getNLogisticPrior(ngrp, pyAnal.prior_cols, pyAnal.restricted)

pyRes = bmdscore.python_nested_result()
pyRes.nparms = pyAnal.parms
pyRes.model = pyAnal.model
bmdsRes = bmdscore.BMDS_results()
boot = bmdscore.nestedBootstrap()
litter = bmdscore.nestedLitterData()
reduced = bmdscore.nestedReducedData()
srData = bmdscore.nestedSRData()
pyRes.bmdsRes = bmdsRes
pyRes.boot = boot
pyRes.litter = litter
pyRes.reduced = reduced
pyRes.srData = srData

bmdscore.pythonBMDSNested(pyAnal, pyRes)


print(pyRes.bmdsRes.BMD)  # noqa: T201
