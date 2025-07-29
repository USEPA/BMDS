import numpy as np
import pytest

from pybmds.datasets.continuous import ContinuousDataset, ContinuousIndividualDataset


class TestFabricateIndividualResponses:
    def setup_method(self):
        self.doses = [0.0, 1.0, 5.0]
        self.ns = [10, 10, 10]
        self.means = [1.0, 2.0, 3.0]
        self.stdevs = [0.1, 0.2, 0.3]

    def test_output_structure(self):
        dataset = ContinuousDataset(
            doses=self.doses, ns=self.ns, means=self.means, stdevs=self.stdevs, name="TestDataset"
        )

        fabricated = dataset.fabricate_individual_responses(seed=42)
        assert isinstance(fabricated, ContinuousIndividualDataset)
        assert len(fabricated.responses) == sum(self.ns)
        assert len(fabricated.individual_doses) == sum(self.ns)

    @pytest.mark.parametrize("tol", [0.001, 0.01, 0.05])
    def test_statistical_accuracy_within_tolerance(self, tol):
        dataset = ContinuousDataset(
            doses=self.doses, ns=self.ns, means=self.means, stdevs=self.stdevs, name="TestDataset"
        )

        fabricated = dataset.fabricate_individual_responses(seed=123, tol=tol)
        for dose, target_mean, target_sd in zip(self.doses, self.means, self.stdevs, strict=False):
            group = np.array(
                [
                    r
                    for d, r in zip(fabricated.individual_doses, fabricated.responses, strict=False)
                    if d == dose
                ]
            )
            calculated_mean = np.mean(group)
            calculated_sd = np.std(group, ddof=0)

            mean_diff = abs(calculated_mean - target_mean)
            sd_diff = abs(calculated_sd - target_sd)

            assert mean_diff < tol, f"Mean mismatch at dose {dose}: {mean_diff} > tol {tol}"
            assert sd_diff < tol, f"SD mismatch at dose {dose}: {sd_diff} > tol {tol}"

    def test_warning_if_cannot_converge(self):
        # Artificially tight tolerance and small max_iter to trigger warning
        dataset = ContinuousDataset(
            doses=[0.0, 1.0, 2.0, 3.0],
            ns=[200, 200, 200, 200],
            means=[1.0, 2.0, 3.0, 4.0],
            stdevs=[2.1, 2.2, 2.3, 2.4],
            name="UnrealisticToleranceTest",
        )
        with pytest.raises(ValueError, match="Could not generate a valid sample"):
            dataset.fabricate_individual_responses(seed=42, tol=1e-10, max_iter=2)

    def test_impose_positivity_enforced(self):
        dataset = ContinuousDataset(
            doses=self.doses, ns=self.ns, means=self.means, stdevs=self.stdevs, name="PositiveOnly"
        )

        fabricated = dataset.fabricate_individual_responses(seed=123, impose_positivity=True)

        assert all(
            r > 0 for r in fabricated.responses
        ), "Found non-positive response despite impose_positivity=True"

    def test_impose_positivity_disabled(self):
        dataset = ContinuousDataset(
            doses=self.doses,
            ns=self.ns,
            means=[0.0, 0.0, 0.0],  # center around 0 to allow negatives
            stdevs=[1.0, 1.0, 1.0],  # wide distribution to increase chance of negatives
            name="AllowNegatives",
        )

        fabricated = dataset.fabricate_individual_responses(seed=456, impose_positivity=False)

        has_negatives = any(r < 0 for r in fabricated.responses)
        assert has_negatives, "Expected some negative values when impose_positivity=False"
