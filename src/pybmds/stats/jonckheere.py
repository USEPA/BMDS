import warnings
from enum import StrEnum
from math import comb

import numpy as np
import scipy.stats as ss
from pydantic import BaseModel

from pybmds.bmdscore import cwilcox

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
        maxVal = 1.0
        minVal = 0.0
        val = -9999
        if self == Hypothesis.two_sided:
            val = 2 * min(decreasing, increasing, 0.5)
        elif self == Hypothesis.increasing:
            val = increasing
        elif self == Hypothesis.decreasing:
            val = decreasing

        if val == -9999:
            raise ValueError("Bad Hypothesis def")  # pragma: no cover

        # ensure p-value is between zero and one
        return min(maxVal, max(minVal, val))


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
                ["Hypothesis", self.hypothesis],
                ["Statistic", self.statistic],
                ["Approach (P-Value)", self.approach],
                ["P-Value", self.p_value],
            ],
            "",
        )


def jonckheere(
    x: np.ndarray,
    group: np.ndarray,
    hypothesis: str | Hypothesis = Hypothesis.two_sided,
    nperm: int | None = None,
    seed: int | None = None,
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
        seed (int, optional): Random number seed used for permutations.

    Returns:
        Test statistic, p-value, alternative hypothesis, and method used.
    """

    # set maximum data length for using exact method
    maxN = 150

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

    statistic = _calc_stat(x, group_indices, group_count)

    if nperm:
        # calculate P-value via permutation
        approach = Approach.permtuation
        rng = np.random.default_rng(seed)
        observed_stat = statistic
        perm_stats = np.empty(nperm, dtype=np.float64)
        x_perm = x.copy()
        for i in range(nperm):
            perm_stats[i] = _calc_stat(x_perm, group_indices, group_count)
            x_perm = rng.permutation(x)
        decreasing = np.sum(perm_stats >= observed_stat) / nperm
        increasing = np.sum(perm_stats <= observed_stat) / nperm

    else:
        # calculate p-value with exact method
        if tot_elements > maxN or len(np.unique(x)) < tot_elements:
            # calculate P-value via normal approximation
            approach = Approach.approximate
            zstat = _zstat(group_size, statistic)
            decreasing = float(ss.norm.cdf(zstat))
            increasing = 1 - decreasing

            if tot_elements < 30:
                msg = "P-value estimated using normal distribution; total observations < 30"
                warnings.warn(msg, stacklevel=2)
        else:
            approach = Approach.exact
            decreasing = sum(_conv_pdf(group_size)[: int(statistic) + 1])
            increasing = 1 - sum(_conv_pdf(group_size)[: int(statistic)])

    pval = hypothesis.calculate_pvalue(decreasing, increasing)

    return TestResult(
        approach=approach, statistic=statistic, p_value=float(pval), hypothesis=hypothesis
    )


def _conv_pdf(group_size: np.ndarray) -> np.ndarray:
    """Finds the probability density function (PDF) using the convolution by Mark van de Wiel.

    Args:
        group_size (np.ndarray): An array of group sizes.

    Returns:
        The probability density function (PDF).
    """
    # Find the cmulative sizes in reverse order and calculate the maximum sum for the PDF computation
    group_count = group_size.size
    csum_groupsize = np.cumsum(np.flip(group_size))
    csum_groupsize = np.flip(csum_groupsize)
    max_sum = sum(group_size[:-1] * csum_groupsize[1:]) + 1

    pdf = np.zeros(max_sum)
    pdf0 = np.zeros(max_sum)
    pdf1 = np.zeros(max_sum)

    # Calculate the PDF for the last two groups
    m = csum_groupsize[-2] - csum_groupsize[-1]
    n = csum_groupsize[-1]
    mn1 = m * n

    for i in range(mn1 + 1):
        pdf[i] = _wilcoxon_pdf(int(i), int(m), int(n))

    # Loop through the remaining groups in reverse order
    for g in range(group_count - 3, -1, -1):
        pdf1[:] = pdf
        pdf.fill(0)

        # Calculate PDF
        m = csum_groupsize[g] - csum_groupsize[g + 1]
        n = csum_groupsize[g + 1]
        mn0 = m * n

        for i in range(mn0 + 1):
            pdf0[i] = _wilcoxon_pdf(int(i), int(m), int(n))

        # Convolve pdf0 and pdf1 and store in pdf
        for i in range(mn0 + 1):
            for j in range(mn1 + 1):
                pdf[i + j] += pdf0[i] * pdf1[j]

        # Update range
        mn1 = mn0 + mn1

    return pdf


def _wilcoxon_pdf(i: int, m: int, n: int) -> np.ndarray:
    """Calculates the raw density for the distribution of the Wilcoxon rank sum statistic.

    Args:
        i (int): Vector of quantiles.
        m: (int): Number of observations in the first sample.
        n: (int): Number of observations in the second sample.

    Returns:
        The raw density for the distribution.
    """
    # Check for NaN values
    if np.isnan(i) or np.isnan(m) or np.isnan(n):
        return i + m + n

    # Check if sample sizes are valid
    if m <= 0 or n <= 0:
        return float("nan")  # Return NaN if invalid

    # Check if i is within valid bounds
    if i < 0 or i > m * n:
        return 0.0  # Return 0 if i is out of bounds

    # Calculate the raw density
    d = cwilcox(i, m, n) / comb(m + n, n)

    return d


def _zstat(group_size: np.ndarray, statistic: int) -> float:
    """Calculation of mean and variance under null hypothesis."""
    ni = group_size[:-1]
    nj = np.cumsum(group_size[::-1])[:-1][::-1]
    j_mean = np.sum(ni * nj / 2)
    j_var = np.sum(ni * nj * (ni + nj + 1) / 12)
    return (statistic - j_mean) / np.sqrt(j_var)


def _calc_stat(x: np.ndarray, group_indices: np.ndarray, group_count: int) -> int:
    """Vectorized Jonckheere statistic calculation."""
    j_sum = 0
    for i in range(group_count - 1):
        xi = x[group_indices == i]
        for j in range(i + 1, group_count):
            xj = x[group_indices == j]
            j_sum += np.sum(xi[:, None] < xj[None, :])
    return j_sum
