# /usr/bin/env python

import numpy as np

# from pybmds import bmdscore
import bmdscore

pyAnal = bmdscore.python_dichotomous_analysis()
pyAnal.model = bmdscore.dich_model.d_weibull
pyAnal.n = 4
pyAnal.Y = np.array([0, 0, 8, 20])
pyAnal.doses = np.array([0, 10, 30, 100])
pyAnal.n_group = np.array([20, 20, 20, 20])
pyAnal.prior = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 1, 1e-6, 18, 18, 100])
pyAnal.BMR = 0.1
pyAnal.BMD_type = 1  # 1 = extra ; added otherwise
pyAnal.alpha = 0.05
pyAnal.degree = 3  # for multistage only
# pyAnal.samples =
# pyAnal.burnin =
pyAnal.parms = 3
pyAnal.prior_cols = 5

pyRes = bmdscore.python_dichotomous_model_result()
pyRes.model = pyAnal.model
pyRes.dist_numE = 200
pyRes.nparms = pyAnal.parms

gof = bmdscore.dichotomous_GOF()
bmdsRes = bmdscore.BMDS_results()
aod = bmdscore.dicho_AOD()
pyRes.gof = gof
pyRes.bmdsRes = bmdsRes
pyRes.aod = aod

bmdscore.pythonBMDSDicho(pyAnal, pyRes)
print(pyRes.bmd)
