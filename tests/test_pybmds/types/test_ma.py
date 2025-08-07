import pybmds


class TestDichotomousModelAverage:
    def test_cpp_str(self, ddataset2):
        # ensure we can generate a string representation of the cpp structs
        session = pybmds.Session(dataset=ddataset2)
        session.add_default_bayesian_models()
        session.execute()
        text = str(session.model_average.structs)
        assert "python_dichotomous_analysis" in text
        assert "python_dichotomousMA_result" in text
        assert len(text.splitlines()) == 2613
