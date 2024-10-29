import numpy as np

import bmdscore


def getMultitumorPrior(degree, prior_cols):
    prG = np.array([0, -17, 0, -18, 18])
    prB = np.array([0, 0.1, 0, 0, 1e4])
    pr = []
    for i in range(prior_cols):
        pr.append(prG[i])
        for _ in range(degree):
            pr.append(prB[i])
    return pr


doses1 = np.array([0, 50, 100, 150, 200])
doses2 = np.array([0, 50, 100, 150, 200])
doses3 = np.array([0, 50, 100, 150, 200])
Y1 = np.array([0, 5, 30, 65, 90])
Y2 = np.array([5, 10, 33, 67, 93])
Y3 = np.array([1, 68, 78, 88, 98])
n_group1 = np.array([100, 100, 100, 100, 100])
n_group2 = np.array([100, 100, 100, 100, 100])
n_group3 = np.array([100, 100, 100, 100, 100])

doses = []
Y = []
n_group = []
doses.append(doses1)
doses.append(doses2)
doses.append(doses3)
Y.append(Y1)
Y.append(Y2)
Y.append(Y3)
n_group.append(n_group1)
n_group.append(n_group2)
n_group.append(n_group3)

dist_numE = 200
ndatasets = 3
BMR = 0.1
BMD_type = 1
alpha = 0.05
prior_cols = 5

n = np.array([5, 5, 5])
degree = np.array([3, 0, 0])
pyAnal = bmdscore.python_multitumor_analysis()
pyAnal.ndatasets = ndatasets
pyAnal.n = n
pyAnal.degree = degree
pyAnal.BMR = BMR
pyAnal.BMD_type = BMD_type
pyAnal.alpha = alpha
pyAnal.prior_cols = prior_cols

pyRes = bmdscore.python_multitumor_result()
pyRes.ndatasets = ndatasets

models = []
resModels = []
nmodels = []
for dataset in range(pyAnal.ndatasets):
    models.append([])
    resModels.append([])
    count = 0
    if degree[dataset] == 0:
        for deg in range(2, pyAnal.n[dataset]):
            models[dataset].append(bmdscore.python_dichotomous_analysis())
            models[dataset][count].model = bmdscore.dich_model.d_multistage
            models[dataset][count].prior = getMultitumorPrior(deg, pyAnal.prior_cols)
            models[dataset][count].degree = deg
            models[dataset][count].parms = deg + 1
            models[dataset][count].Y = Y[dataset]
            models[dataset][count].n_group = n_group[dataset]
            models[dataset][count].doses = doses[dataset]
            models[dataset][count].n = len(Y[dataset])
            models[dataset][count].BMR = BMR
            models[dataset][count].BMD_type = BMD_type
            models[dataset][count].alpha = alpha
            models[dataset][count].prior_cols = prior_cols
            resModels[dataset].append(bmdscore.python_dichotomous_model_result())
            resModels[dataset][count].model = bmdscore.dich_model.d_multistage
            resModels[dataset][count].nparms = deg + 1
            resModels[dataset][count].dist_numE = dist_numE
            gof = bmdscore.dichotomous_GOF()
            bmdsRes = bmdscore.BMDS_results()
            aod = bmdscore.dicho_AOD()
            resModels[dataset][count].gof = gof
            resModels[dataset][count].bmdsRes = bmdsRes
            resModels[dataset][count].aod = aod
            count = count + 1
    else:
        models[dataset].append(bmdscore.python_dichotomous_analysis())
        models[dataset][count].model = bmdscore.dich_model.d_multistage
        models[dataset][count].prior = getMultitumorPrior(degree[dataset], pyAnal.prior_cols)
        models[dataset][count].degree = degree[dataset]
        models[dataset][count].parms = degree[dataset] + 1
        models[dataset][count].Y = Y[dataset]
        models[dataset][count].n_group = n_group[dataset]
        models[dataset][count].doses = doses[dataset]
        models[dataset][count].n = len(Y[dataset])
        models[dataset][count].BMR = BMR
        models[dataset][count].BMD_type = BMD_type
        models[dataset][count].alpha = alpha
        models[dataset][count].prior_cols = prior_cols
        resModels[dataset].append(bmdscore.python_dichotomous_model_result())
        resModels[dataset][count].model = bmdscore.dich_model.d_multistage
        resModels[dataset][count].nparms = degree[dataset] + 1
        resModels[dataset][count].dist_numE = dist_numE
        gof = bmdscore.dichotomous_GOF()
        bmdsRes = bmdscore.BMDS_results()
        aod = bmdscore.dicho_AOD()
        resModels[dataset][count].gof = gof
        resModels[dataset][count].bmdsRes = bmdsRes
        resModels[dataset][count].aod = aod
        count = 1
    nmodels.append(count)


pyAnal.models = models
pyRes.models = resModels
pyAnal.nmodels = nmodels
pyRes.nmodels = nmodels

bmdscore.pythonBMDSMultitumor(pyAnal, pyRes)
print(pyRes.BMD)  # noqa: T201
