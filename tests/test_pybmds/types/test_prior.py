from copy import deepcopy
from textwrap import dedent

import numpy as np
import pytest

from pybmds.constants import DistType, PriorClass, PriorDistribution
from pybmds.models.continuous import Polynomial
from pybmds.models.dichotomous import Multistage
from pybmds.types.priors import ModelPriors, Prior


@pytest.fixture
def mock_prior():
    t = PriorDistribution.Uniform
    return ModelPriors(
        prior_class=PriorClass.frequentist_restricted,
        priors=[
            Prior(name="a", type=t, initial_value=1, stdev=1, min_value=1, max_value=1),
            Prior(name="b", type=t, initial_value=2, stdev=2, min_value=2, max_value=2),
            Prior(name="c", type=t, initial_value=3, stdev=3, min_value=3, max_value=3),
        ],
        variance_priors=[
            Prior(name="d", type=t, initial_value=4, stdev=4, min_value=4, max_value=4),
            Prior(name="e", type=t, initial_value=5, stdev=5, min_value=5, max_value=5),
        ],
    )


@pytest.fixture
def mock_nested_dichotomous_prior():
    t = PriorDistribution.Uniform
    return ModelPriors(
        prior_class=PriorClass.frequentist_restricted,
        priors=[
            Prior(name="a", type=t, initial_value=2, stdev=0, min_value=1, max_value=4),
            Prior(name="b", type=t, initial_value=3, stdev=0, min_value=2, max_value=5),
            Prior(name="phi", type=t, initial_value=4, stdev=0, min_value=3, max_value=6),
        ],
        variance_priors=[],
    )


class TestModelPriors:
    def test_get_prior(self, mock_prior):
        assert mock_prior.get_prior("a").name == "a"
        assert mock_prior.get_prior("d").name == "d"
        with pytest.raises(ValueError):
            mock_prior.get_prior("z")

    def test_update(self, mock_prior):
        a = mock_prior.get_prior("a")
        assert a.initial_value == 1
        mock_prior.update("a", initial_value=2)
        assert a.initial_value == 2

    def test_to_c(self, mock_prior):
        # fmt: off
        assert np.allclose(
            mock_prior.to_c(),
            [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
        )
        assert np.allclose(
            mock_prior.to_c(degree=1),
            [0.0, 0.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0]
        )
        assert np.allclose(
            mock_prior.to_c(degree=2),
            [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
        )
        assert np.allclose(
            mock_prior.to_c(degree=3),
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 1.0, 2.0, 3.0, 3.0, 1.0, 2.0, 3.0, 3.0, 1.0, 2.0, 3.0, 3.0]
        )
        assert np.allclose(
            mock_prior.to_c(degree=4),
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 1.0, 2.0, 3.0, 3.0, 3.0, 1.0, 2.0, 3.0, 3.0, 3.0, 1.0, 2.0, 3.0, 3.0, 3.0]
        )
        assert np.allclose(mock_prior.to_c(
            dist_type=DistType.normal),
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 5.0, 1.0, 2.0, 3.0, 5.0, 1.0, 2.0, 3.0, 5.0, 1.0, 2.0, 3.0, 5.0]
        )
        assert np.allclose(
            mock_prior.to_c(dist_type=DistType.log_normal),
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 5.0, 1.0, 2.0, 3.0, 5.0, 1.0, 2.0, 3.0, 5.0, 1.0, 2.0, 3.0, 5.0]
        )
        assert np.allclose(
            mock_prior.to_c(dist_type=DistType.normal_ncv),
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        )
        # fmt: on

    def test_to_c_nd(self, mock_nested_dichotomous_prior):
        prior = mock_nested_dichotomous_prior
        assert np.allclose(prior.to_c_nd(n_phi=1), [1, 2, 3, 4, 5, 6])
        assert np.allclose(prior.to_c_nd(n_phi=2), [1, 2, 3, 3, 4, 5, 6, 6])
        assert np.allclose(prior.to_c_nd(n_phi=3), [1, 2, 3, 3, 3, 4, 5, 6, 6, 6])

    def test_multistage_update(self, ddataset):
        m = Multistage(dataset=ddataset, settings=dict(degree=8))
        expected = dedent(
            """
            ╒═════════════╤═══════════╤═══════╤═══════╕
            │ Parameter   │   Initial │   Min │   Max │
            ╞═════════════╪═══════════╪═══════╪═══════╡
            │ g           │         0 │   -18 │    18 │
            │ b1          │         0 │     0 │ 10000 │
            │ b2          │         0 │     0 │ 10000 │
            │ b3          │         0 │     0 │ 10000 │
            │ b4          │         0 │     0 │ 10000 │
            │ b5          │         0 │     0 │ 10000 │
            │ b6          │         0 │     0 │ 10000 │
            │ b7          │         0 │     0 │ 10000 │
            │ b8          │         0 │     0 │ 10000 │
            ╘═════════════╧═══════════╧═══════╧═══════╛
        """
        )
        assert m.priors_tbl() == expected.strip()
        m = Multistage(dataset=ddataset, settings=dict(degree=8))
        m.settings.priors.update("g", min_value=-0.1, initial_value=0.05, max_value=0.1)
        m.settings.priors.update("b1", max_value=10)
        m.settings.priors.update("b2", max_value=2)
        m.settings.priors.update("b3", max_value=3)
        m.settings.priors.update("b4", max_value=4)
        m.settings.priors.update("b5", max_value=5)
        m.settings.priors.update("b6", max_value=6)
        m.settings.priors.update("b7", max_value=7)
        m.settings.priors.update("b8", max_value=8)
        expected = dedent(
            """
            ╒═════════════╤═══════════╤═══════╤═══════╕
            │ Parameter   │   Initial │   Min │   Max │
            ╞═════════════╪═══════════╪═══════╪═══════╡
            │ g           │      0.05 │  -0.1 │   0.1 │
            │ b1          │      0    │   0   │  10   │
            │ b2          │      0    │   0   │   2   │
            │ b3          │      0    │   0   │   3   │
            │ b4          │      0    │   0   │   4   │
            │ b5          │      0    │   0   │   5   │
            │ b6          │      0    │   0   │   6   │
            │ b7          │      0    │   0   │   7   │
            │ b8          │      0    │   0   │   8   │
            ╘═════════════╧═══════════╧═══════╧═══════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

    def test_polynomial_update(self, cdataset):
        m = Polynomial(dataset=cdataset, settings=dict(degree=8))
        expected = dedent(
            """
            ╒═════════════╤═══════════╤═════════╤═══════╕
            │ Parameter   │   Initial │     Min │   Max │
            ╞═════════════╪═══════════╪═════════╪═══════╡
            │ g           │         0 │   0     │  1000 │
            │ b1          │         0 │ -18     │     0 │
            │ b2          │         0 │  -1e+06 │     0 │
            │ b3          │         0 │  -1e+06 │     0 │
            │ b4          │         0 │  -1e+06 │     0 │
            │ b5          │         0 │  -1e+06 │     0 │
            │ b6          │         0 │  -1e+06 │     0 │
            │ b7          │         0 │  -1e+06 │     0 │
            │ b8          │         0 │  -1e+06 │     0 │
            │ rho         │         0 │   0     │    18 │
            │ alpha       │         0 │ -18     │    18 │
            ╘═════════════╧═══════════╧═════════╧═══════╛
        """
        )
        assert m.priors_tbl() == expected.strip()
        m = Polynomial(dataset=cdataset, settings=dict(degree=8))
        m.settings.priors.update("g", min_value=-0.1, initial_value=0.05, max_value=0.1)
        m.settings.priors.update("b1", min_value=-1, max_value=1)
        m.settings.priors.update("b2", max_value=2)
        m.settings.priors.update("b3", max_value=3)
        m.settings.priors.update("b4", max_value=4)
        m.settings.priors.update("b5", max_value=5)
        m.settings.priors.update("b6", max_value=6)
        m.settings.priors.update("b7", max_value=7)
        m.settings.priors.update("b8", max_value=8)
        expected = dedent(
            """
            ╒═════════════╤═══════════╤═════════╤═══════╕
            │ Parameter   │   Initial │     Min │   Max │
            ╞═════════════╪═══════════╪═════════╪═══════╡
            │ g           │      0.05 │  -0.1   │   0.1 │
            │ b1          │      0    │  -1     │   1   │
            │ b2          │      0    │  -1e+06 │   2   │
            │ b3          │      0    │  -1e+06 │   3   │
            │ b4          │      0    │  -1e+06 │   4   │
            │ b5          │      0    │  -1e+06 │   5   │
            │ b6          │      0    │  -1e+06 │   6   │
            │ b7          │      0    │  -1e+06 │   7   │
            │ b8          │      0    │  -1e+06 │   8   │
            │ rho         │      0    │   0     │  18   │
            │ alpha       │      0    │ -18     │  18   │
            ╘═════════════╧═══════════╧═════════╧═══════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # assert that full JSON validation and parsing works as expected
        initial = m.settings.priors.model_dump_json()
        assert ModelPriors.model_validate_json(initial).model_dump_json() == initial

    def test_nested_dichotomous_update(self, mock_nested_dichotomous_prior):
        prior = mock_nested_dichotomous_prior
        prior.update("phi1", min_value=-10, max_value=10)
        assert np.allclose(prior.to_c_nd(n_phi=2), [1, 2, -10, 3, 4, 5, 10, 6])

    def test_priors_list_validation(self, mock_prior):
        # valid if fixed
        priors = deepcopy(mock_prior)
        priors.update("a", min_value=1, initial_value=1, max_value=1)
        assert len(priors.priors_list()) > 0

        for settings, message in [
            (dict(min_value=1, max_value=0), "Min Value > Max Value"),
            (dict(min_value=1, initial_value=0), "Initial Value < Min Value"),
            (dict(max_value=1, initial_value=2), "Initial Value > Max Value"),
        ]:
            priors = deepcopy(mock_prior)
            priors.update("a", **settings)
            with pytest.warns(UserWarning, match=message):
                priors.priors_list()
