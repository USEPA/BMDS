from math import comb

import numpy as np
from pydantic import BaseModel
from scipy.stats import norm

from ..constants import ZEROISH
from ..utils import pretty_table


class TestResult(BaseModel):
    statistic: float
    p_value_asymptotic: float
    p_value_exact: float

    def tbl(self) -> str:
        return pretty_table(
            [
                ["Statistic", self.statistic],
                ["P-Value (Asymptotic)", self.p_value_asymptotic],
                ["P-Value (Exact)", self.p_value_exact],
            ],
            "",
        )


def cochran_armitage(dose: np.ndarray, n: np.ndarray, incidence: np.ndarray) -> TestResult:
    """
    Cochran-Armitage trend test for response data across ordered dose levels.

    This function computes both the asymptotic and conditional exact p-values for testing
    the presence of a monotonic trend in incidence rates across increasing dose groups.

    The asymptotic p-value is based on a normal approximation of the linear trend statistic
    proposed by Cochran (1954) and Armitage (1955). The the exact one-sided p-value for the
    Cochran-Armitage trend test using a special case of the linear rank test algorithm by
    Mehta, Patel, and Tsiatis (1992).

    Args:
        dose (np.ndarray): independent variable, monotonically increasing
        n (np.ndarray): number per group
        incidence (np.ndarray): number affected per group

    References:
        Cochran, W. G. (1954). Some methods for strengthening the common χ² tests.
        Biometrics, 10(4), 417-451.

        Armitage, P. (1955). Tests for linear trends in proportions and frequencies.
        Biometrics, 11(3), 375-386.

        Mehta, C. R., Patel, N. R., & Tsiatis, A. A. (1992). Exact stratified linear rank tests
        for ordered categorical and binary data. Journal of Computational and Graphical
        Statistics, 1(1), 21-40.
    """
    # Input validation
    if dose.size < 3:
        raise ValueError("At least three dose groups are required.")
    if not (dose.size == n.size == incidence.size):
        raise ValueError("All input vectors must be of the same length.")
    if np.any(n <= 1):
        raise ValueError("All dose groups must have at least two individuals.")
    if np.any(incidence < 0):
        raise ValueError("The number of responders must be nonnegative.")
    if np.any(incidence > n):
        raise ValueError("The number of responders cannot exceed group total.")
    if not np.all(np.diff(dose) > 0):
        raise ValueError("Dose levels must be strictly increasing.")

    total_n = n.sum()
    total_inc = incidence.sum()
    pbar = total_inc / total_n
    weights = n / total_n

    numerator = np.sum((incidence - pbar * n) * dose)
    denominator = np.sqrt(np.sum(weights * dose**2) - (np.sum(weights * dose)) ** 2)
    factor = np.sqrt(total_n / ((total_n - total_inc) * total_inc))
    z_statistic = -factor * numerator / denominator

    asymptotic_pvalue = norm.cdf(z_statistic)
    exact_pvalue = _p_value_exact(dose, n, incidence)

    return TestResult(
        statistic=z_statistic,
        p_value_asymptotic=float(asymptotic_pvalue),
        p_value_exact=exact_pvalue,
    )


def _p_value_exact(dose: np.ndarray, n: np.ndarray, incidence: np.ndarray) -> float:
    """
    Compute the exact one-sided p-value for the Cochran-Armitage trend test
    using a special case of the linear rank test algorithm by Mehta, Patel, and Tsiatis (1992).

    This implementation computes the conditional exact distribution of a linear trend statistic
    for a single 2xc contingency table with fixed marginal totals, corresponding to a test for
    monotonic trend in binary response data across ordered dose groups.

    The algorithm uses a network-based recursive traversal of the response space,
    enumerating all possible response configurations that match the observed marginals,
    and computes the proportion of those configurations with test statistics
    as extreme or more extreme than the observed value.

    References:
        Mehta, C. R., Patel, N. R., & Tsiatis, A. A. (1992). Exact stratified linear rank tests
        for ordered categorical and binary data. Journal of Computational and Graphical
        Statistics, 1(1), 21-40.
    """
    dk = dose / dose[1]  # normalize
    rest = dk - np.floor(dk)
    rest = rest[rest > 0]
    if len(rest) > 0:
        mult = min(1 / np.prod(rest), 1e12)
    else:
        mult = 1
    dk = np.round(dk * mult).astype(int)

    tot_k = len(dk)
    total_n = np.sum(n)
    total_inc = np.sum(incidence)
    target_score = np.sum(incidence * dk)

    nodes = [[] for _ in range(tot_k + 1)]
    nodes[0] = [0]
    for i in range(tot_k):
        low = max(0, total_inc - np.sum(n[i + 1 :]))
        high = min(total_inc, np.sum(n[: i + 1]))
        nodes[i + 1] = list(range(low, high + 1))

    # Build arcs: (from, to, score_contrib, count_paths)
    arcs = [[] for _ in range(tot_k)]
    for i in range(tot_k):
        for j in nodes[i]:
            for k in nodes[i + 1]:
                if k >= j and k - j <= n[i]:
                    score = dk[i] * (k - j)
                    count = comb(n[i], k - j)
                    arcs[i].append((j, k, score, count))

    score_paths = {0: {0: 1}}

    for i in range(tot_k):
        next_score_paths = {}
        for j, k, score, count in arcs[i]:
            if j not in score_paths:
                continue
            for s, c in score_paths[j].items():
                new_s = s + score
                new_c = c * count
                if k not in next_score_paths:
                    next_score_paths[k] = {}
                if new_s in next_score_paths[k]:
                    next_score_paths[k][new_s] += new_c
                else:
                    next_score_paths[k][new_s] = new_c
        score_paths = next_score_paths

    total_comb = comb(total_n, total_inc)

    passed = sum(
        c for s, c in score_paths.get(total_inc, {}).items() if s >= target_score - ZEROISH
    )

    return passed / total_comb
