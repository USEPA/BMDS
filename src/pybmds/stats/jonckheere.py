import warnings
from enum import StrEnum

import numpy as np
import scipy.stats as ss
from pydantic import BaseModel

from ..utils import pretty_table


class Hypothesis(StrEnum):
    two_sided = "two-sided"
    increasing = "increasing"
    decreasing = "decreasing"

    def calculate_pvalue(self, decreasing: float, increasing: float) -> float:
        """Compute the p-value with the given alternative hypothesis.

        Args:
            decreasing (float): The probability of observing a value less than a certain threshold.
            increasing (float): The probability of observing a value greater than a certain threshold.
        Returns:
            P-value with given alternative hypothesis.
        """
        if self == Hypothesis.two_sided:
            return 2 * min(decreasing, increasing, 0.5)
        elif self == Hypothesis.increasing:
            return increasing
        elif self == Hypothesis.decreasing:
            return decreasing
        raise ValueError("Unreachable code path")  # pragma: no cover


class Approach(StrEnum):
    exact = "exact"
    approximate = "approximate"
    permtuation = "permutation"


class TestResult(BaseModel):
    approach: Approach
    statistic: float
    p_value: float
    hypothesis: Hypothesis

    def tbl(self) -> str:
        return pretty_table(
            [
                ["Approach", self.approach],
                ["Hypothesis", self.hypothesis],
                ["Statistic", self.statistic],
                ["P Value", self.p_value],
            ],
            "",
        )


def jonckheere(
    x: np.ndarray,
    group: np.ndarray,
    hypothesis: str | Hypothesis = Hypothesis.two_sided,
    nperm: int | None = None,
) -> TestResult:
    """Compute Jonckheere-Terpstra test.

    This function computes the Jonckheere-Terpstra trend test which is a non-parametric statistical
    test used to detect a trend in data across multiple ordered groups. The test statistic counts
    the number of times values in higher order groups are greater than those in lower order
    groups. When not using permutations, the p-value is estimated using a normal approximation.

    Args:
        x (np.ndarray): A numeric array of data values.
        group (np.ndarray): A numeric or ordered factor indicating group membership.
        hypothesis (str | Hypothesis): The type of hypothesis test to be run (two-sided, increasing, or decreasing).
        nperm (int, optional): The number of permutations for the test.

    Returns:
        Test statistic, p-value, alternative hypothesis, and method used.
    """
    # Input validation
    if not np.issubdtype(x.dtype, np.number):
        raise ValueError("Data needs to be numeric")
    if not (np.issubdtype(group.dtype, np.number) or np.issubdtype(group.dtype, object)):
        raise ValueError("Group needs to be numeric or ordered factor")
    if len(group) != len(x):
        raise ValueError("Data and group values need to be the same length")

    # will raise ValueError if not valid choice
    hypothesis = Hypothesis(hypothesis)

    # Remove missing data using numpy vectorization
    mask = np.isfinite(x) & np.isfinite(group)
    x = x[mask]
    group = group[mask]

    tot_elements = x.size
    if tot_elements == 0:
        raise ValueError("Either data or group is missing for all observations")

    # Use numpy for group sizes and sorting
    unique_groups, group_indices = np.unique(group, return_inverse=True)
    group_count = unique_groups.size
    group_size = np.bincount(group_indices)

    if group_count <= 1:
        raise ValueError("Only one group has non-missing data")

    # Vectorized calculation of Jonckheere statistic
    # For each pair of groups (i < j), count number of times x_i < x_j
    j_sum = 0
    for i in range(group_count - 1):
        xi = x[group_indices == i]
        for j in range(i + 1, group_count):
            xj = x[group_indices == j]
            j_sum += np.sum(xi[:, None] < xj[None, :])

    statistic = j_sum

    if nperm:
        # calculate P-value via permutation
        pval = _perm(x, group_indices, group_size, hypothesis, nperm)
        approach = Approach.permtuation

    else:
        # calculate P-value via normal approximation
        approach = Approach.approximate
        zstat = _zstat(group_size, statistic)
        decreasing = float(ss.norm.cdf(zstat))
        increasing = 1 - decreasing
        pval = hypothesis.calculate_pvalue(decreasing, increasing)

        if tot_elements < 30:
            msg = "P-value estimated using normal distribution; total observations < 30"
            warnings.warn(msg, stacklevel=2)

    return TestResult(
        approach=approach, statistic=statistic, p_value=float(pval), hypothesis=hypothesis
    )


def _zstat(group_size: np.ndarray, statistic: int) -> float:
    """Calculation of mean and variance under null hypothesis."""
    ni = group_size[:-1]
    nj = np.cumsum(group_size[::-1])[:-1][::-1]
    j_mean = np.sum(ni * nj / 2)
    j_var = np.sum(ni * nj * (ni + nj + 1) / 12)
    return (statistic - j_mean) / np.sqrt(j_var)


def _perm(
    x: np.ndarray,
    group_indices: np.ndarray,
    group_size: np.ndarray,
    hypothesis: Hypothesis,
    nperm: int,
) -> float:
    """Fast permutation test using numpy vectorization."""
    rng = np.random.default_rng()
    observed_stat = _calc_stat(x, group_indices, group_size)
    perm_stats = np.empty(nperm, dtype=np.float64)
    for i in range(nperm):
        perm_stats[i] = _calc_stat(x, group_indices, group_size)
        x = rng.permutation(x)
    decreasing = np.sum(perm_stats >= observed_stat) / nperm
    increasing = np.sum(perm_stats <= observed_stat) / nperm
    return hypothesis.calculate_pvalue(decreasing, increasing)


def _calc_stat(
    x: np.ndarray,
    group_indices: np.ndarray,
    group_size: np.ndarray,
) -> int:
    """Vectorized Jonckheere statistic calculation."""
    csum_groupsize = np.concatenate(([0], np.cumsum(group_size)))
    group_count = len(group_size)
    stat = 0
    for i in range(group_count - 1):
        xi = x[csum_groupsize[i] : len(x)]
        rank_a = ss.rankdata(xi)
        rank_b = rank_a[: group_size[i]]
        stat += sum(rank_b) - group_size[i] * (group_size[i] + 1) / 2
    return stat
