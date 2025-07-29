import numpy as np

_ERR_MSG = "Could not generate a valid sample for target_mean and target_sd."


def exact_rnorm(
    n: int,
    target_mean: float,
    target_sd: float,
    seed: int | None = 42,
    tolerance: float = 0.01,
    max_iterations: int = 10_000,
    impose_positivity: bool = True,
) -> tuple[np.ndarray, int]:
    """Calculate a sample of values within the tolerance of the distribution summary statistics.

    Args:
        n (int): Number of samples to generate.
        target_mean (float): The target mean of the samples.
        target_sd (float): The target standard deviation fo the samples.
        seed (int, optional): Random number seed; defaults to 42.
        tolerance (float, optional): Maximum difference between simulated and target values.
        max_iterations (int, optional): _Maximum iterations to attempt.
        impose_positivity (bool, optional): Require values to be positive.

    Raises:
        ValueError: _description_

    Returns:
        tuple[np.ndarray, int]: _description_
    """
    if n == 1:
        if target_sd != 0:
            raise ValueError(_ERR_MSG)
        return np.array([target_mean]), 0
    if n == 2:
        if impose_positivity and target_mean - target_sd < 0:
            raise ValueError(_ERR_MSG)
        return np.array([target_mean - target_sd, target_mean + target_sd]), 0

    rng = np.random.default_rng(seed)
    closest_diff = float("inf")
    for iteration in range(1, max_iterations + 1):
        x = rng.normal(size=n)
        x = (x - np.mean(x)) / np.std(x) * target_sd + target_mean

        if impose_positivity and x[x >= 0].size < n:
            continue

        mean_diff = abs(np.mean(x) - target_mean)
        sd_diff = abs(np.std(x) - target_sd)
        current_diff = mean_diff + sd_diff

        if current_diff < closest_diff:
            closest_diff = current_diff

        if (
            mean_diff / abs(target_mean + 1e-8) < tolerance
            and sd_diff / abs(target_sd + 1e-8) < tolerance
        ):
            return x, iteration

    raise ValueError(_ERR_MSG)
