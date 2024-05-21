from math import log
from typing import Self

from pydantic import BaseModel
from scipy import stats


class Test(BaseModel):
    DF: float = -1
    CDF: float = -1
    SS: float = 0.0
    MSE: float = 0.0
    AIC: float = 0.0
    TEST: float = -1.0


class AnovaTests(BaseModel):
    test1: Test
    test2: Test
    test3: Test

    @staticmethod
    def compute_likelihoods(n_obs, ns, ym, yd):
        Ntot = ns[0]
        Ntot = sum(ns)
        sigma2 = yd[0] * (ns[0] - 1)
        for i in range(1, n_obs):
            sigma2 += yd[i] * (ns[i] - 1)
        sigma2 = sigma2 / Ntot
        lkA1 = -Ntot * (1.0 + log(sigma2)) / 2.0
        lkA2 = -ns[0] * log(yd[0] * (ns[0] - 1) / ns[0]) / 2.0 - Ntot / 2.0
        for i in range(1, n_obs):
            lkA2 -= ns[i] * log(yd[i] * (ns[i] - 1) / ns[i]) / 2.0
        ybar = ym[0] * ns[0]
        for i in range(1, n_obs):
            ybar += ym[i] * ns[i]
        ybar = ybar / Ntot
        sigma2 = yd[0] * (ns[0] - 1) + ns[0] * (ym[0] - ybar) * (ym[0] - ybar)
        for i in range(1, n_obs):
            sigma2 += yd[i] * (ns[i] - 1) + ns[i] * (ym[i] - ybar) * (ym[i] - ybar)
        sigma2 = sigma2 / Ntot
        lkR = -Ntot * (1.0 + log(sigma2)) / 2.0
        lkA3 = lkA1
        return lkA1, lkA2, lkA3, lkR

    @classmethod
    def get_anova_c3_tests(cls, nparm, n_obs, a1, a2, a3, ar) -> Self:
        # Based on Hill model, modified from DTMS3ANOVAC.c

        # The xlk is not real
        xlk = (a1 + a2 + a3 + ar) / 4
        parm_known = 1
        anovas = [Test() for i in range(5)]

        # Compute DF and assign LLK for each test
        anovas[0].DF = n_obs + 1
        anovas[0].SS = a1

        anovas[1].DF = 2 * n_obs
        anovas[1].SS = a2

        anovas[2].DF = n_obs + 2 - parm_known
        anovas[2].SS = a3

        anovas[3].DF = 2
        anovas[3].SS = ar

        anovas[4].DF = nparm - 2
        anovas[4].SS = xlk

        # Compute likelihood ratio MSE and CDF
        anovas[0].MSE = 2 * (a2 - a1)
        anovas[0].CDF = anovas[1].DF - anovas[0].DF

        anovas[1].MSE = 2 * (a2 - a3)
        anovas[1].CDF = anovas[1].DF - anovas[2].DF

        anovas[2].MSE = 2 * (a3 - xlk)
        anovas[2].CDF = anovas[2].DF - anovas[4].DF

        anovas[3].MSE = 2 * (a2 - ar)
        anovas[3].CDF = anovas[1].DF - anovas[3].DF

        for anova in anovas:
            anova.AIC = -2 * (anova.SS - anova.DF)
            if anova.MSE >= 0.0 and anova.CDF > 0:
                anova.TEST = 1 - stats.chi2.cdf(anova.MSE, anova.CDF)

        # Only return test 1, test 2 and test 3 in the order
        return cls(test1=anovas[3], test2=anovas[0], test3=anovas[1])

    @staticmethod
    def output_3tests(tests) -> str:
        if tests is None:
            return "ANOVA cannot be calculated for this dataset."

        outputs = [
            "                     Tests of Interest    ",
            "   Test    -2*log(Likelihood Ratio)  Test df        P-Value    ",
        ]
        for i, test in enumerate([tests.test1, tests.test2, tests.test3]):
            outputs.append("   Test %d %20.6g %10d %16.4g" % (i + 1, test.MSE, test.CDF, test.TEST))
        return "\n".join(outputs)
