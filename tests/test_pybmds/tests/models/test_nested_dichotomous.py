from textwrap import dedent

from pybmds.constants import PriorClass
from pybmds.models import nested_dichotomous


class TestNestedLogistic:
    def test_param_names(self, nd_dataset, nd_dataset4):
        assert nd_dataset.num_dose_groups == 3
        assert nd_dataset4.num_dose_groups == 4

        model = nested_dichotomous.NestedLogistic(nd_dataset)
        names = model.get_param_names()
        assert names == ["g", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3"]

        model = nested_dichotomous.NestedLogistic(nd_dataset4)
        names = model.get_param_names()
        assert names == ["g", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3", "phi4"]

    def test_model_priors(self, nd_dataset):
        m = nested_dichotomous.NestedLogistic(dataset=nd_dataset)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤════════╤════════╕
        │ Parameter   │   Initial │    Min │    Max │
        ╞═════════════╪═══════════╪════════╪════════╡
        │ g           │         0 │  0     │  1     │
        │ b           │         0 │ -1e+06 │  1e+06 │
        │ theta1      │         0 │  0     │  1     │
        │ theta2      │         0 │ -1e+06 │  1e+06 │
        │ rho         │         0 │  1     │ 18     │
        │ phi1        │         0 │  0     │  1e+06 │
        │ phi2        │         0 │  0     │  1e+06 │
        │ phi3        │         0 │  0     │  1e+06 │
        ╘═════════════╧═══════════╧════════╧════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check overrides
        m.settings.priors.update("g", min_value=1, max_value=2)
        m.settings.priors.update("phi2", min_value=5, max_value=6)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤════════╤════════╕
        │ Parameter   │   Initial │    Min │    Max │
        ╞═════════════╪═══════════╪════════╪════════╡
        │ g           │         0 │  1     │  2     │
        │ b           │         0 │ -1e+06 │  1e+06 │
        │ theta1      │         0 │  0     │  1     │
        │ theta2      │         0 │ -1e+06 │  1e+06 │
        │ rho         │         0 │  1     │ 18     │
        │ phi1        │         0 │  0     │  1e+06 │
        │ phi2        │         0 │  5     │  6     │
        │ phi3        │         0 │  0     │  1e+06 │
        ╘═════════════╧═══════════╧════════╧════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check prior classes
        m = nested_dichotomous.NestedLogistic(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_restricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 1
        m = nested_dichotomous.NestedLogistic(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_unrestricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 0

    def test_execute(self, nd_dataset4):
        # add seed for reproducibility
        analysis = nested_dichotomous.NestedLogistic(nd_dataset4, settings=dict(bootstrap_seed=1))
        analysis.execute()
        text = analysis.text()
        assert len(text) > 0


class TestNctr:
    def test_param_names(self, nd_dataset, nd_dataset4):
        assert nd_dataset.num_dose_groups == 3
        assert nd_dataset4.num_dose_groups == 4

        model = nested_dichotomous.Nctr(nd_dataset)
        names = model.get_param_names()
        assert names == ["g", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3"]

        model = nested_dichotomous.Nctr(nd_dataset4)
        names = model.get_param_names()
        assert names == ["g", "b", "theta1", "theta2", "rho", "phi1", "phi2", "phi3", "phi4"]

    def test_model_priors(self, nd_dataset):
        m = nested_dichotomous.Nctr(dataset=nd_dataset)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤═════════╤════════════╕
        │ Parameter   │   Initial │     Min │        Max │
        ╞═════════════╪═══════════╪═════════╪════════════╡
        │ g           │         0 │  0      │  18        │
        │ b           │         0 │  0      │ -18        │
        │ theta1      │         0 │ -0.0625 │  -0.142857 │
        │ theta2      │         0 │ -0.0625 │  -0.142857 │
        │ rho         │         0 │  1      │  18        │
        │ phi1        │         0 │  0      │  18        │
        │ phi2        │         0 │  0      │  18        │
        │ phi3        │         0 │  0      │  18        │
        ╘═════════════╧═══════════╧═════════╧════════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check overrides
        m.settings.priors.update("g", min_value=1, max_value=2)
        m.settings.priors.update("phi2", min_value=5, max_value=6)
        expected = dedent(
            """
        ╒═════════════╤═══════════╤═════════╤════════════╕
        │ Parameter   │   Initial │     Min │        Max │
        ╞═════════════╪═══════════╪═════════╪════════════╡
        │ g           │         0 │  1      │   2        │
        │ b           │         0 │  0      │ -18        │
        │ theta1      │         0 │ -0.0625 │  -0.142857 │
        │ theta2      │         0 │ -0.0625 │  -0.142857 │
        │ rho         │         0 │  1      │  18        │
        │ phi1        │         0 │  0      │  18        │
        │ phi2        │         0 │  5      │   6        │
        │ phi3        │         0 │  0      │  18        │
        ╘═════════════╧═══════════╧═════════╧════════════╛
        """
        )
        assert m.priors_tbl() == expected.strip()

        # check prior classes
        m = nested_dichotomous.Nctr(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_restricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 1
        m = nested_dichotomous.Nctr(
            dataset=nd_dataset, settings=dict(priors=PriorClass.frequentist_unrestricted)
        )
        assert m.settings.priors.get_prior("rho").min_value == 0
