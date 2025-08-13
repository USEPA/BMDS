from enum import StrEnum

import numpy as np
import pandas as pd
import scipy.stats as ss
from pydantic import BaseModel
from scipy.special import comb


class Hypothesis(StrEnum):
    two_sided = "two-sided"
    increasing = "increasing"
    decreasing = "decreasing"

    def calc_pval(self, decreasing: float, increasing: float) -> float:
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


class TestResult(BaseModel):
    statistic: float
    p_value: float
    hypothesis: Hypothesis


def jonckheere(
    x: np.ndarray,
    group: np.ndarray,
    hypothesis: str | Hypothesis = Hypothesis.two_sided,
    nperm: int | None = None,
) -> TestResult:
    """Compute Jonckheere-Terpstra test.

    This function computes the Jonckheere-Terpstra trend test
    which is a non-parametric statistical test used to detect
    a trend in data across multiple ordered groups.

    Args:
        x (np.ndarray): A numeric array of data values.
        group (np.ndarray): A numeric or ordered factor indicating group membership.
        hypothesis (str | Hypothesis): The type of hypothesis test to be run (two-sided, increasing, or decreasing).
        nperm (int, optional): The numer of permutations for the test.

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

    x = x[np.isfinite(x) & np.isfinite(group)]

    tot_elements = x.size
    if tot_elements == 0:
        raise ValueError("Either data or group is missing for all observations")

    group_size = pd.Series(group).value_counts()
    group_count = group_size.size
    group_size = group_size.to_numpy()

    if group_count <= 1:
        raise ValueError("Only one group has non-missing data")

    # Find cumulative sum of the group sizes
    csum_groupsize = np.concatenate(([0], np.cumsum(group_size)))

    j_sum = 0
    j_mean = 0
    j_var = 0

    # Loop through each group to find the test statistics
    for i in range(group_count - 1):
        current_size = group_size[i]
        remain_size = tot_elements - csum_groupsize[i + 1]
        j_sum += (
            np.sum(ss.rankdata(x[csum_groupsize[i] : tot_elements])[:current_size])
            - current_size * (current_size + 1) / 2
        )
        j_mean += current_size * remain_size / 2
        j_var += current_size * remain_size * (current_size + remain_size + 1) / 12

    # Adjust for increasing/decreasing trends
    j_sum = 2 * j_mean - j_sum
    statistic = j_sum

    # Calculate p-value based on permutation or normal approximation
    j_index = int(j_sum)
    if nperm:
        pval = _jtperm(x, group_count, group_size, csum_groupsize, hypothesis, nperm)
    else:
        if tot_elements > 100 or len(np.unique(x)) < tot_elements:
            zstat = (statistic - j_mean) / np.sqrt(j_var)
            decreasing = float(ss.norm.cdf(zstat))
            increasing = 1 - decreasing
            pval = hypothesis.calc_pval(decreasing, increasing)
        else:
            decreasing = sum(_conv_pdf(group_size)[: j_index + 1])
            increasing = 1 - sum(_conv_pdf(group_size)[:j_index])
            pval = hypothesis.calc_pval(decreasing, increasing)

    return TestResult(statistic=statistic, p_value=float(pval), hypothesis=hypothesis)


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


def _jtperm(
    x: np.ndarray,
    group_count: int,
    group_size: np.ndarray,
    csum_groupsize: np.ndarray,
    hypothesis: Hypothesis,
    nperm: int,
) -> float:
    """Calculates the p-value with the given number of permuations.

    Args:
        x (np.ndarray): A numeric array of data values.
        group_count (int): The number of groups.
        group_size (np.ndarray): Array indicating the size of each group.
        csum_groupsize (np.ndarray): Array containing the cumulative sum of the group sizes.
        hypothesis (str | Hypothesis): The type of hypothesis test to be run (two-sided, increasing, or decreasing).
        nperm (int): The numer of permutations for the test.

    Returns:
        The p-value given the number of permutations specified.
    """
    perm_jsum = np.zeros(nperm)
    rng = np.random.default_rng()
    x_copy = x.copy()
    # Calculate the test statistic for each group for the specified number of permutations then store results
    for j in range(nperm):
        j_sum = 0
        for i in range(group_count - 1):
            current_size = group_size[i]
            ranks = np.argsort(np.argsort(x_copy[csum_groupsize[i] :]))
            j_sum += np.sum(ranks[:current_size]) - current_size * (current_size + 1) / 2
        perm_jsum[j] = j_sum
        rng.shuffle(x_copy)

    decreasing = np.sum(perm_jsum >= perm_jsum[0]) / nperm
    increasing = np.sum(perm_jsum <= perm_jsum[0]) / nperm

    return hypothesis.calc_pval(decreasing, increasing)


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
    d = _cwilcox(i, m, n) / comb(m + n, n)

    return d


def _cwilcox(k, m, n):
    """Calculates the number of ways to obtain a sum from a set of numbers given the constraints of the Wilcoxon rank sum test.

    Args:
        k (int): Sum to achieve
        m: (int): Max value for each element
        n: (int): Max number of elements

    Returns:
        The number of ways to obtain the sum with the given inputs.
    """
    # Empty vectors for storing results
    w = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    # Calculate the total number of observations
    u = m * n

    # Check if k is too large or too small
    if k < 0 or k > u:
        return 0

    # Calculate half of the total observations
    c = u // 2

    # Adjust k if it's greater than half
    if k > c:
        k = u - k

    # Ensure m is the smaller sample size
    i, j = (m, n) if m < n else (n, m)

    # Check if j is 0
    if j == 0:
        return k == 0

    # If k is less than the size of the second sample, recurse
    if j > 0 and k < j:
        return _cwilcox(k, i, k)

    # Initialize the array if not already done
    if w[i - 1][j - 1] is None:
        w[i - 1][j - 1] = [-1] * (c + 1)

    # If the value has not been computed yet
    if w[i - 1][j - 1][k - 1] < 0:
        if j == 0:
            w[i - 1][j - 1][k - 1] = k == 0
        else:
            # Recursive calls to compute the value
            w[i - 1][j - 1][k - 1] = _cwilcox(k - j, i - 1, j) + _cwilcox(k, i, j - 1)

    return w[i - 1][j - 1][k - 1]
