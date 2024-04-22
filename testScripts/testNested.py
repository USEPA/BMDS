# /usr/bin/env python

# from pybmds import bmdscore
import bmdscore
import numpy as np

pyAnal = bmdscore.python_nested_analysis()
pyAnal.model = bmdscore.nested_model.nlogistic
pyAnal.restricted = True
pyAnal.doses = np.array([])
pyAnal.litterSize = np.array([])
pyAnal.incidence = np.array([])
pyAnal.lsc = np.array([])
pyAnal.LSC_type = 1
pyAnal.ILC_type = 1
pyAnal.BMD_type = 1
pyAnal.background = 1
pyAnal.BMR = 0.1
pyAnal.alpha = 0.05
pyAnal.iterations = 1000
pyAnal.seed = -9999

pyRes = bmdscore.python_nested_result()
bmdsRes = bmdscore.BMDS_results()
boot = bmdscore.nestedBootstrap()
litter = bmdscore.nestedLitterData()
reduced = bmdscore.nestedReducedData()
pyRes.bmdsRes = bmdsRes
pyRes.boot = boot
pyRes.litter = litter
pyRes.reduced = reduced

bmdscore.pythonBMDSNested(pyAnal, pyRes)


print(pyRes.bmdsRes.BMD)  # noqa: T201
