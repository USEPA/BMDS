from enum import StrEnum

import numpy as np
import pandas as pd
import scipy.stats as ss
from pydantic import BaseModel
from scipy.stats import norm


class Alternative(StrEnum):
    two_sided = "two-sided"
    increasing = "increasing"
    decreasing = "decreasing"

    def calc_pval(self, decreasing: float, increasing: float) -> float:
        if self == Alternative.two_sided:
            return 2 * min(decreasing, increasing, 0.5)
        elif self == Alternative.increasing:
            return increasing
        elif self == Alternative.decreasing:
            return decreasing
        raise ValueError("Unreachable code path")  # pragma: no cover


class JonckheereResult(BaseModel):
    statistic: float
    p_value: float
    alternative: Alternative


def jonckheere(
    x, group, alternative: str | Alternative = Alternative.two_sided, nperm: int | None = None
) -> JonckheereResult:
    """
    Compute Jonckheere-Terpstra test.

    Returns test statistic, p-value, alternative hypothesis, and method used.
    """

    if not np.issubdtype(x.dtype, np.number):
        raise ValueError("Data needs to be numeric")
    if not (np.issubdtype(group.dtype, np.number) or np.issubdtype(group.dtype, object)):
        raise ValueError("Group needs to be numeric or ordered factor")
    if len(group) != len(x):
        raise ValueError("Data and group values need to be the same length")

    # will raise ValueError if not valid choice
    alternative = Alternative(alternative)

    x = x[np.isfinite(x) & np.isfinite(group)]

    n = x.size
    if n == 0:
        raise ValueError("Either data or group is missing for all observations")

    group_size = pd.Series(group).value_counts()
    group_count = group_size.size
    group_size = group_size.to_numpy()

    if group_count <= 1:
        raise ValueError("Only one group has non-missing data")

    csum_groupsize = np.concatenate(([0], np.cumsum(group_size)))

    jtrsum = 0
    jtmean = 0
    jtvar = 0

    for i in range(group_count - 1):
        na = group_size[i]
        nb = n - csum_groupsize[i + 1]
        jtrsum += np.sum(ss.rankdata(x[csum_groupsize[i] : n])[:na]) - na * (na + 1) / 2
        jtmean += na * nb / 2
        jtvar += na * nb * (na + nb + 1) / 12

    jtrsum = int(2 * jtmean - jtrsum)
    statistic = jtrsum

    if nperm:
        pval = _jtperm(x, group_count, group_size, csum_groupsize, alternative, nperm)
    else:
        if n > 100 or len(np.unique(x)) < n:
            zstat = (statistic - jtmean) / np.sqrt(jtvar)
            decreasing = ss.norm.cdf(zstat)
            increasing = 1 - decreasing
            pval = alternative.calc_pval(decreasing, increasing)
        else:
            decreasing = sum(_conv_pdf(group_size)[1 : (jtrsum + 1)])
            increasing = 1 - sum(_conv_pdf(group_size)[1:(jtrsum)])
            pval = alternative.calc_pval(decreasing, increasing)

    return JonckheereResult(statistic=statistic, p_value=pval, alternative=alternative)


def _conv_pdf(group_size: np.ndarray) -> np.ndarray:
    """
    Find the probability density function (PDF) using the convolution by Mark van de Wiel.
    """
    group_count = group_size.size
    csum_groupsize = np.cumsum(np.flip(group_size))
    csum_groupsize = np.flip(csum_groupsize)
    max_sum = sum(group_size[:-1] * csum_groupsize[1:]) + 1
    return _jon_pdf(max_sum, group_count, csum_groupsize)


def _jtperm(
    x: np.ndarray,
    group_count: int,
    group_size: np.ndarray,
    csum_groupsize,
    alternative: Alternative,
    nperm: int,
) -> float:
    """
    Evaluate p-value.
    """

    pjtrsum = np.zeros(nperm)

    for j in range(nperm):
        jtrsum = 0
        for i in range(group_count - 1):
            na = group_size[i]
            ranks = np.argsort(np.argsort(x[csum_groupsize[i] :]))
            jtrsum += np.sum(ranks[:na]) - na * (na + 1) / 2

        pjtrsum[j] = jtrsum

    decreasing = np.sum(pjtrsum >= pjtrsum[0]) / nperm
    increasing = np.sum(pjtrsum <= pjtrsum[0]) / nperm
    return alternative.calc_pval(decreasing, increasing)


def _wilcox_normal_pdf(x: float, m: float, n: float) -> np.ndarray:
    mu = m * (m + n + 1) / 2
    sigma = np.sqrt(m * n * (m + n + 1) / 12)
    return norm.pdf(x, loc=mu, scale=sigma)


def _jon_pdf(max_sum: int, group_count: int, csum_groupsize: np.ndarray) -> np.ndarray:
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
