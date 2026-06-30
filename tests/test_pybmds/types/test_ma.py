from types import SimpleNamespace

import numpy as np

import pybmds
from pybmds.constants import BMDS_BLANK_VALUE
from pybmds.types.ma import DichotomousModelAverageResult


class TestDichotomousModelAverage:
    def test_cpp_str(self, ddataset2):
        # ensure we can generate a string representation of the cpp structs
        session = pybmds.Session(dataset=ddataset2)
        session.add_default_bayesian_models()
        session.execute()
        text = str(session.model_average.structs)
        assert "python_dichotomous_analysis" in text
        assert "python_dichotomousMA_result" in text


class TestDichotomousModelAverageResultHelpers:
    def test_loud_weights_and_posterior_fallbacks(self):
        result = DichotomousModelAverageResult
        assert np.all(result._loud_weights([np.nan, BMDS_BLANK_VALUE]) == BMDS_BLANK_VALUE)
        assert result._loud_weights([1.0, 2.0]).sum() == 1.0
        assert result._loud_posteriors(np.array([-1.0, -1.0]), 1, [1.0, 2.0], []).sum() == 1.0
        averaged = result._loud_posteriors(np.array([-1.0, -1.0]), 0, [1.0, 2.0], [2.0, 1.0])
        assert averaged.tolist() == [0.5, 0.5]

    def test_draw_shape_helpers(self):
        result = DichotomousModelAverageResult
        draws = np.arange(6.0)
        assert result._reshape_flat_draws(draws, 4).shape == (6,)
        assert result._reshape_flat_draws(draws, 2).shape == (2, 3)

        loud = SimpleNamespace(BMD=[1.0, 2.0], parms=[1.0, 2.0, 3.0, 4.0])
        bmd, parms = result._loud_draws(loud, 1)
        assert bmd.shape == (2,)
        assert parms.shape == (2, 2)

    def test_combined_loud_result_falls_back_to_per_chain_results(self):
        chains = [SimpleNamespace(BMD=[1.0])]
        cpp_result = SimpleNamespace(combinedLoudRes=None, loudRes=chains)
        assert DichotomousModelAverageResult._combined_loud_result(cpp_result) is chains
