import numpy as np

import pybmds
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
            bmd_dist=np.array([[1.0, np.nan], [2.0, 3.0]]),
            priors=np.array([1.0]),
            posteriors=np.array([1.0]),
            model_bmd_dist=[np.array([[1.0, np.inf], [2.0, 3.0]])],
            model_parm_dist=[np.array([[[1.0], [np.nan]], [[2.0], [3.0]]])],
            dr_x=np.array([0.0, 1.0]),
            dr_y=np.array([0.0, 0.5]),
        )

        data = result.model_dump()

        assert data["bmd_dist"] == [[1.0, None], [2.0, 3.0]]
        assert data["model_bmd_dist"] == [[[1.0, None], [2.0, 3.0]]]
        assert data["model_parm_dist"] == [[[[1.0], [None]], [[2.0], [3.0]]]]

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
