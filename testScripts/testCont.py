# /usr/bin/env python

import bmdscore
import numpy as np

# from pybmds import bmdscore

pyAnal = bmdscore.python_continuous_analysis()
pyAnal.model = bmdscore.cont_model.exp_5
# isIncreasing needed if detectAdvDir is set to false
# pyAnal.isIncreasing
pyAnal.n = 5
pyAnal.doses = np.array([0, 62.5, 125, 250, 500])
pyAnal.Y = np.array([24.3, 27, 31.4, 39.3, 54.2])
pyAnal.n_group = np.array([10, 10, 10, 10, 10])
pyAnal.sd = np.array([4.93, 3.16, 7.05, 13.2, 25.8])
pyAnal.BMD_type = 2
pyAnal.BMR = 1.0
pyAnal.suff_stat = True
pyAnal.parms = 5
pyAnal.prior_cols = 5
pyAnal.transform_dose = 0
pyAnal.prior = np.array(
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1, 0.5, 0.2, 1, 0, 0, -20, 1, -18, 100, 100, 20, 18, 18]
)
pyAnal.disttype = bmdscore.distribution.normal
pyAnal.alpha = 0.05
pyAnal.detectAdvDir = True
pyAnal.restricted = True

pyRes = bmdscore.python_continuous_model_result()
pyRes.model = pyAnal.model
pyRes.dist_numE = 200
pyRes.nparms = pyAnal.parms

gof = bmdscore.continuous_GOF()
bmdsRes = bmdscore.BMDS_results()
aod = bmdscore.continuous_AOD()
aod.TOI = bmdscore.testsOfInterest()
pyRes.gof = gof
pyRes.bmdsRes = bmdsRes
pyRes.aod = aod

bmdscore.pythonBMDSCont(pyAnal, pyRes)
print(pyRes.bmd)  # noqa: T201
