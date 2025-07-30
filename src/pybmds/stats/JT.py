import numpy as np
import pandas as pd
import scipy.stats as ss
from scipy.stats import norm


def jonckheere(x, group, alternative="two.sided", nperm=None):
    """
    This function runs the Jonckheere-Terpstra test and returns
    the test statistic, p-value, alternative hypothesis, and method
    used
    """

    if not np.issubdtype(x.dtype, np.number):
        raise ValueError("Data needs to be numeric")
    if not (np.issubdtype(group.dtype, np.number) or np.issubdtype(group.dtype, object)):
        raise ValueError("Group needs to be numeric or ordered factor")
    if len(group) != len(x):
        raise ValueError("Data and group values need to be the same length")

    alternative_opt = ["two.sided", "increasing", "decreasing"]
    if alternative not in alternative_opt:
        raise ValueError("Alternative choice not valid")

    finite_filter = np.isfinite(x) & np.isfinite(group)
    x = x[finite_filter]
    n = len(x)

    if n == 0:
        raise ValueError("Either data or group is missing for all observations")

    group_size = pd.Series(group).value_counts()
    group_count = len(group_size)
    group_size = group_size.to_numpy()

    if group_count <= 1:
        raise ValueError("Only one group has non-missing data")

    csum_groupsize = np.concatenate(([0], np.cumsum(group_size)))

    jtrsum = jtmean = jtvar = 0

    for i in range(group_count - 1):
        na = group_size[i]
        nb = n - csum_groupsize[i + 1]
        jtrsum += np.sum(ss.rankdata(x[csum_groupsize[i] : n])[:na]) - na * (na + 1) / 2
        jtmean += na * nb / 2
        jtvar += na * nb * (na + nb + 1) / 12

    jtrsum = int(2 * jtmean - jtrsum)
    STATISTIC = jtrsum

    if nperm is not None:
        PVAL = _jtperm(x, group_count, group_size, csum_groupsize, alternative, nperm)
        pass
    else:
        if n > 100 or len(np.unique(x)) < n:
            zstat = (STATISTIC - jtmean) / np.sqrt(jtvar)
            PVAL = ss.norm.cdf(zstat)
            if alternative == "two.sided":
                PVAL = 2 * min(PVAL, 1 - PVAL, 0.5)
            elif alternative == "increasing":
                PVAL = 1 - PVAL
            elif alternative == "decreasing":
                PVAL = PVAL
        else:
            dPVAL = sum(_conv_pdf(group_size)[1 : (jtrsum + 1)])
            iPVAL = 1 - sum(_conv_pdf(group_size)[1:(jtrsum)])
            if alternative == "two.sided":
                PVAL = 2 * min(iPVAL, dPVAL, 0.5)
            elif alternative == "increasing":
                PVAL = iPVAL
            elif alternative == "decreasing":
                PVAL = dPVAL

            pass

    return {
        "statistic": STATISTIC,
        "p.value": PVAL,
        "alternative": alternative,
        "method": "Jonckheere-Terpstra test",
    }


def _conv_pdf(group_size):
    """
    This function finds the probability density
    function (PDF) using the convolution by
    Mark van de Wiel
    """

    group_count = len(group_size)
    csum_groupsize = np.cumsum(np.flip(group_size))
    csum_groupsize = np.flip(csum_groupsize)
    max_sum = sum(group_size[:-1] * csum_groupsize[1:]) + 1
    pdf = _jon_pdf(max_sum, group_count, csum_groupsize)
    return pdf


def _jtperm(x, group_count, group_size, csum_groupsize, alternative, nperm):
    """
    This function evalvulates p-value
    """

    pjtrsum = np.zeros(nperm)

    for j in range(nperm):
        jtrsum = 0
        for i in range(group_count - 1):
            na = group_size[i]
            ranks = np.argsort(np.argsort(x[csum_groupsize[i] :]))
            jtrsum += np.sum(ranks[:na]) - na * (na + 1) / 2

        pjtrsum[j] = jtrsum
        np.random.Generator(x)

    iPVAL = np.sum(pjtrsum <= pjtrsum[0]) / nperm
    dPVAL = np.sum(pjtrsum >= pjtrsum[0]) / nperm

    if alternative == "two.sided":
        PVAL = 2 * min(iPVAL, dPVAL, 0.5)
    elif alternative == "increasing":
        PVAL = iPVAL
    elif alternative == "decreasing":
        PVAL = dPVAL

    return PVAL


def _wilcox_normal_pdf(x, m, n):
    mu = m * (m + n + 1) / 2
    sigma = np.sqrt(m * n * (m + n + 1) / 12)
    return norm.pdf(x, loc=mu, scale=sigma)


def _jon_pdf(max_sum, group_count, csum_groupsize):
    pdf = np.zeros(max_sum)
    pdf0 = np.zeros(max_sum)
    pdf1 = np.zeros(max_sum)
    group_count_index = group_count - 1

    m = csum_groupsize[group_count_index - 1] - csum_groupsize[group_count_index]
    n = csum_groupsize[group_count_index]
    mn1 = m * n

    for i in range(mn1 + 1):
        pdf[i] = _wilcox_normal_pdf(float(i), float(m), float(n))

    for g in range(group_count_index - 2, 0, -1):
        pdf1[:] = pdf
        pdf.fill(0)

        m = csum_groupsize[g] - csum_groupsize[g + 1]
        n = csum_groupsize[g + 1]
        mn0 = m * n

        for i in range(mn0):
            pdf0[i] = _wilcox_normal_pdf(float(i), float(m), float(n))

        for i in range(mn0):
            for j in range(mn0):
                pdf[i + j] += pdf0[i] * pdf1[j]

        mn1 = mn0 + mn1

    return pdf
