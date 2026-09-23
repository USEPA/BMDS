from types import SimpleNamespace

import numpy as np

import pybmds
from pybmds.constants import BMDS_BLANK_VALUE
from pybmds.types.cma import ContinuousModelAverageResult


class TestContinuousModelAverage:
    def test_cpp_str(self, cdataset3):
        # ensure we can generate a string representation of the cpp structs
        session = pybmds.Session(dataset=cdataset3)
        session.add_default_bayesian_models()
        session.execute()
        text = str(session.model_average.structs)
        assert "python_continuous_analysis" in text
        assert "python_continuousMA_result" in text

    def test_loud_json_draws_preserve_shape_with_nulls(self):
        result = ContinuousModelAverageResult(
            bmdl=1.0,
            bmd=2.0,
            bmdu=3.0,
            bmdl_y=0.1,
            bmd_y=0.2,
            bmdu_y=0.3,
            bmd_dist=np.array([[1.0, np.nan], [BMDS_BLANK_VALUE, 3.0]]),
            priors=np.array([1.0]),
            posteriors=np.array([1.0]),
            model_bmd_dist=[np.array([[1.0, np.inf], [BMDS_BLANK_VALUE, 3.0]])],
            model_parm_dist=[np.array([[[1.0], [np.nan]], [[BMDS_BLANK_VALUE], [3.0]]])],
            dr_x=np.array([0.0, 1.0]),
            dr_y=np.array([0.0, 0.5]),
        )

        data = result.model_dump()

        assert data["bmd_dist"] == [[1.0, None], [None, 3.0]]
        assert data["model_bmd_dist"] == [[[1.0, None], [None, 3.0]]]
        assert data["model_parm_dist"] == [[[[1.0], [None]], [[None], [3.0]]]]

    def test_without_loud_draws_removes_raw_draw_arrays(self):
        result = ContinuousModelAverageResult(
            bmdl=1.0,
            bmd=2.0,
            bmdu=3.0,
            bmdl_y=0.1,
            bmd_y=0.2,
            bmdu_y=0.3,
            bmd_dist=np.array([[1.0, 2.0], [3.0, 4.0]]),
            priors=np.array([1.0]),
            posteriors=np.array([1.0]),
            model_bmd_dist=[np.array([[1.0, 2.0]])],
            model_parm_dist=[np.array([[[1.0], [2.0]]])],
            dr_x=np.array([0.0, 1.0]),
            dr_y=np.array([0.0, 0.5]),
        )

        trimmed = result.without_loud_draws()

        assert trimmed.model_dump()["bmd_dist"] == []
        assert trimmed.model_dump()["model_bmd_dist"] == []
        assert trimmed.model_dump()["model_parm_dist"] == []
        assert result.model_dump()["bmd_dist"] == [[1.0, 2.0], [3.0, 4.0]]

    def test_draw_shape_helpers_accept_per_chain_transposed_parameters(self):
        chains = [
            SimpleNamespace(
                BMD=[1.0, 2.0, 3.0],
                parms=np.array([[10.0, 20.0, 30.0], [11.0, 21.0, 31.0]]),
            ),
            SimpleNamespace(
                BMD=[4.0, 5.0, 6.0],
                parms=np.array([[40.0, 50.0, 60.0], [41.0, 51.0, 61.0]]),
            ),
        ]

        bmd, parms = ContinuousModelAverageResult._loud_draws(chains, 2)

        assert bmd.shape == (2, 3)
        assert parms.shape == (2, 3, 2)
        np.testing.assert_array_equal(parms[0], [[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]])

    def test_combined_loud_result_prefers_combined_draws_and_falls_back(self):
        combined = SimpleNamespace(BMD=[1.0], parms=[[2.0]])
        chains = [SimpleNamespace(BMD=[3.0], parms=[[4.0]])]

        assert (
            ContinuousModelAverageResult._combined_loud_result(
                SimpleNamespace(combinedLoudRes=combined, loudRes=chains)
            )
            is combined
        )
        assert (
            ContinuousModelAverageResult._combined_loud_result(
                SimpleNamespace(combinedLoudRes=SimpleNamespace(BMD=[]), loudRes=chains)
            )
            is chains
        )
