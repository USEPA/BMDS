import json

import pytest

import pybmds
from pybmds.constants import Models, PriorClass


class TestSession:
    def test_dichotomous(self, ddataset2, rewrite_data_files, data_path):
        # make sure serialize looks correct
        session1 = pybmds.Session(id=1, name="test", description="hello", dataset=ddataset2)
        session1.add_default_models()
        session1.execute_and_recommend()
        d = session1.to_dict()

        # spot check a few keys
        # -> dataset
        assert d["dataset"]["doses"] == [0.0, 50.0, 100.0, 150.0, 200.0]
        # -> models (with results)
        assert len(d["models"]) == 10
        assert list(d["models"][0].keys()) == ["name", "model_class", "settings", "results"]
        # -> models average
        assert d["model_average"] is None
        # -> models recommendation
        assert d["recommender"]["settings"]["enabled"] is True
        assert d["recommender"]["results"]["recommended_model_variable"] == "aic"
        assert d["selected"]["model_index"] is None

        # ensure we can convert back to a session from JSON serialization
        session2 = pybmds.Session.from_serialized(json.loads(json.dumps(d)))
        assert isinstance(session2, pybmds.Session)
        assert session2.dataset.doses == [0.0, 50.0, 100.0, 150.0, 200.0]
        assert len(session2.models) == 10
        assert session2.models[0].has_results is True

        # make sure we get the same result back after deserializing
        d1 = session1.serialize().model_dump()
        d2 = session2.serialize().model_dump()
        assert d1 == d2

        # dataframe
        df = session1.to_df()
        assert "slope_factor" not in df.columns

        # docx
        docx = session1.to_docx(session_inputs_table=True)

        if rewrite_data_files:
            df.to_excel(data_path / "reports/session-dichotomous.xlsx", index=False)
            docx.save(data_path / "reports/session-dichotomous.docx")

    def test_dichotomous_ma(self, ddataset2, data_path, rewrite_data_files):
        # make sure serialize looks correct
        session1 = pybmds.Session(dataset=ddataset2)
        session1.add_default_models()
        session1.add_model_averaging()
        session1.execute_and_recommend()
        d = session1.to_dict()

        if rewrite_data_files:
            (data_path / "reports/session-dichotomous.json").write_text(
                session1.serialize().model_dump_json()
            )

        # spot check a few keys
        assert d["dataset"]["doses"] == [0.0, 50.0, 100.0, 150.0, 200.0]
        assert len(d["models"]) == 10
        assert list(d["models"][0].keys()) == ["name", "model_class", "settings", "results"]
        assert d["model_average"]["model_indexes"] == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert "bmd" in d["model_average"]["results"]

        # ensure we can convert back to a session
        session2 = pybmds.Session.from_serialized(json.loads(json.dumps(d)))
        assert isinstance(session2, pybmds.Session)
        assert session2.dataset.doses == [0.0, 50.0, 100.0, 150.0, 200.0]
        assert len(session2.models) == 10
        assert session2.models[0].has_results is True

        # make sure we get the same result back after deserializing
        d1 = session1.serialize().model_dump()
        d2 = session2.serialize().model_dump()
        assert d1 == d2

    @pytest.mark.skip(reason="TODO - fix pybmds nested dichotomous binding")
    def test_nested_dichotomous(self, nd_dataset4, rewrite_data_files, data_path):
        session = pybmds.Session(dataset=nd_dataset4)
        session.add_default_models()
        session.execute_and_recommend()

        # check serialization
        d = session.to_dict()
        session2 = session.from_serialized(d)
        assert session.to_dict() == session2.to_dict()

        # dataframe
        df = session.to_df()

        # docx
        docx = session.to_docx(session_inputs_table=True, all_models=True)

        if rewrite_data_files:
            df.to_excel(data_path / "reports/session-nested-dichotomous.xlsx", index=False)
            docx.save(data_path / "reports/session-nested-dichotomous.docx")

    def test_dll_version(self, ddataset2):
        session = pybmds.Session(dataset=ddataset2)
        version = session.dll_version()
        assert isinstance(version, str)
        assert int(version.split(".")[0]) >= 2021  # assume dll in format "YYYY.MM..."


class TestSessionPlot:
    @pytest.mark.mpl_image_compare
    def test_dichotomous_colorize(self, ddataset):
        session = pybmds.Session(dataset=ddataset)
        session.add_model(Models.Weibull)
        session.add_model(Models.Logistic)
        session.execute()
        return session.plot(colorize=True)

    @pytest.mark.mpl_image_compare
    def test_dichotomous(self, ddataset):
        session = pybmds.Session(dataset=ddataset)
        session.add_model(Models.Weibull)
        session.add_model(Models.Logistic)
        session.execute()
        return session.plot(colorize=False)

    @pytest.mark.mpl_image_compare
    def test_dichotomous_bayesian_colorize(self, ddataset):
        session = pybmds.Session(dataset=ddataset)
        session.add_model(Models.Weibull, {"priors": PriorClass.bayesian})
        session.add_model(Models.Logistic, {"priors": PriorClass.bayesian})
        session.execute()
        return session.plot(colorize=True)

    @pytest.mark.mpl_image_compare
    def test_dichotomous_bayesian(self, ddataset):
        session = pybmds.Session(dataset=ddataset)
        session.add_model(Models.Weibull, {"priors": PriorClass.bayesian})
        session.add_model(Models.Logistic, {"priors": PriorClass.bayesian})
        session.execute()
        return session.plot(colorize=False)

    @pytest.mark.mpl_image_compare
    def test_dichotomous_ma_colorize(self, ddataset):
        session = pybmds.Session(dataset=ddataset)
        session.add_model(Models.Weibull, {"priors": PriorClass.bayesian})
        session.add_model(Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model_averaging()
        session.execute()
        return session.plot(colorize=True)

    @pytest.mark.mpl_image_compare
    def test_dichotomous_ma(self, ddataset):
        session = pybmds.Session(dataset=ddataset)
        session.add_model(Models.Weibull, {"priors": PriorClass.bayesian})
        session.add_model(Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model_averaging()
        session.execute()
        return session.plot(colorize=False)

    @pytest.mark.mpl_image_compare
    def test_continuous_colorize(self, cdataset):
        session = pybmds.Session(dataset=cdataset)
        session.add_model(Models.Linear)
        session.add_model(Models.ExponentialM3)
        session.execute()
        return session.plot(colorize=True)

    @pytest.mark.mpl_image_compare
    def test_continuous(self, cdataset):
        session = pybmds.Session(dataset=cdataset)
        session.add_model(Models.Linear)
        session.add_model(Models.ExponentialM3)
        session.execute()
        return session.plot(colorize=False)

    @pytest.mark.mpl_image_compare
    def test_continuous_individual_colorize(self, cidataset):
        session = pybmds.Session(dataset=cidataset)
        session.add_model(Models.Linear)
        session.add_model(Models.ExponentialM3)
        session.execute()
        return session.plot(colorize=True)

    @pytest.mark.mpl_image_compare
    def test_continuous_individual(self, cidataset):
        session = pybmds.Session(dataset=cidataset)
        session.add_model(Models.Linear)
        session.add_model(Models.ExponentialM3)
        session.execute()
        return session.plot(colorize=False)
