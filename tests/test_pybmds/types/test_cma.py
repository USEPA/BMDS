import pybmds


class TestContinuousModelAverage:
    def test_cpp_str(self, cdataset2):
        # ensure we can generate a string representation of the cpp structs
        session = pybmds.Session(dataset=cdataset2)
        session.add_default_bayesian_models()
        session.execute()
        text = str(session.model_average.structs)
        assert "python_continuous_analysis" in text
        assert "python_continuousMA_result" in text
