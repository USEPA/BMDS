#/usr/bin/env python

import numpy as np
#from pybmds import bmdscore
import bmdscore

pyMA = bmdscore.python_dichotomousMA_analysis()
pyMA.nmodels = 9
pyMA.actual_parms=np.array([4,3,2,3,3,4,2,2,3])
pyMA.prior_cols = np.full((pyMA.nmodels),5)
pyMA.models = np.array([1,2,3,4,5,6,7,8,9])
pyMA.modelPriors = np.full(pyMA.nmodels, 1.0/pyMA.nmodels)
priorList = []
p1=np.array([1, 1, 1, 2, -1, 0, -3, 0.693147, 2, 3, 3.3, 0.5, -40, -40, -40, 0, 40, 40, 40, 40])
p2=np.array([1, 2, 2, 0, 0.693147, 0, 2, 0.424264, 1, -18, 0.2, 0, 18, 20, 10000])
p3=np.array([1,2,0, 0, 2, 2, -20, 0, 20, 40])
p4=np.array([1,1,2,0,0,0.693147,2,1,0.5,-20,-40,0,20,40,20])
p5=np.array([1,1,2,0,0,0.693147,2,1,0.5,-20,-40,0,20,40,20])
p6=np.array([1,2,2,2,0,0,0,0,2,1,1,1,-20,0,0,0,20,1e6,1e6,1e6])
p7=np.array([1,2,0,0,2,1,-20,0,20,40])
p8=np.array([1,2,0,0,2,1,-20,0,20,18])
p9=np.array([1,2,2,0,0.424264,0,2,0.5,1.5,-20,0,0,20,40,10000])


priorList.append(p1)
priorList.append(p2)
priorList.append(p3)
priorList.append(p4)
priorList.append(p5)
priorList.append(p6)
priorList.append(p7)
priorList.append(p8)
priorList.append(p9)
pyMA.priors = priorList
#pyMA.priors=np.asarray(priorList)

pyAnal = bmdscore.python_dichotomous_analysis()
pyAnal.BMD_type = 1
pyAnal.BMR = 0.1
pyAnal.alpha = 0.05
pyAnal.Y= np.array([2, 10, 13, 15, 15])
pyAnal.n_group = np.array([14, 15, 15, 15, 15])
pyAnal.doses = np.array([0, 11, 30, 100, 356])
pyAnal.n = 5
pyMA.pyDA = pyAnal

pyMARes = bmdscore.python_dichotomousMA_result()
pyMARes.nmodels = pyMA.nmodels
pyMARes.dist_numE = 200

resList = []
r1 = bmdscore.python_dichotomous_model_result() 
r2 = bmdscore.python_dichotomous_model_result() 
r3 = bmdscore.python_dichotomous_model_result() 
r4 = bmdscore.python_dichotomous_model_result() 
r5 = bmdscore.python_dichotomous_model_result()  
r6 = bmdscore.python_dichotomous_model_result() 
r7 = bmdscore.python_dichotomous_model_result() 
r8 = bmdscore.python_dichotomous_model_result() 
r9 = bmdscore.python_dichotomous_model_result() 

resList.append(r1)
resList.append(r2)
resList.append(r3)
resList.append(r4)
resList.append(r5)
resList.append(r6)
resList.append(r7)
resList.append(r8)
resList.append(r9)
pyMARes.models = resList

for x in range(pyMA.nmodels):
   pyMARes.models[x].model = pyMA.models[x]
   pyMARes.models[x].nparms = pyMA.actual_parms[x]
   pyMARes.models[x].dist_numE = pyMARes.dist_numE

bmdsRes = bmdscore.BMDSMA_results()
bmdsRes.BMD = np.full(pyMA.nmodels, -9999)
bmdsRes.BMDL = np.full(pyMA.nmodels, -9999)
bmdsRes.BMDU = np.full(pyMA.nmodels, -9999)
bmdsRes.ebUpper = np.full(pyAnal.n, -9999)
bmdsRes.ebLower = np.full(pyAnal.n, -9999)

pyMARes.bmdsRes = bmdsRes

bmdscore.pythonBMDSDichoMA(pyMA, pyAnal, pyMARes, bmdsRes)



 
