import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pybmds
from pybmds.constants import DistType, Models, PriorClass


class TestSession:
    def test_add_default_models_continuous_uses_legacy_defaults_for_non_loud_priors(
        self, cdataset3
    ):
        session = pybmds.Session(dataset=cdataset3)

        session.add_default_models({"priors": PriorClass.frequentist_restricted})

        model_types = [type(model) for model in session.models]
        assert len(session.models) == 7
        assert pybmds.models.continuous.Linear in model_types
        assert model_types.count(pybmds.models.continuous.Polynomial) == 2
        assert pybmds.models.continuous.MultiplicativeHill not in model_types
        assert pybmds.models.continuous.InverseExponential not in model_types
        assert pybmds.models.continuous.Lognormal not in model_types
        assert pybmds.models.continuous.Gamma not in model_types
        assert pybmds.models.continuous.LMS not in model_types

    def test_add_default_models_continuous_excludes_loud_only_models(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)

        session.add_default_models({"priors": PriorClass.bayesian_loud})

        model_types = {type(model) for model in session.models}
        assert len(session.models) == 4
        assert pybmds.models.continuous.Linear not in model_types
        assert pybmds.models.continuous.Polynomial not in model_types
        assert pybmds.models.continuous.MultiplicativeHill not in model_types
        assert pybmds.models.continuous.InverseExponential not in model_types
        assert pybmds.models.continuous.Lognormal not in model_types
        assert pybmds.models.continuous.Gamma not in model_types
        assert pybmds.models.continuous.LMS not in model_types

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
        docx = session1.to_docx(session_inputs_table=True, all_models=True, bmd_cdf_table=True)

        if rewrite_data_files:
            df.to_excel(data_path / "reports/session-dichotomous.xlsx", index=False)
            docx.save(data_path / "reports/session-dichotomous.docx")

    def test_dichotomous_ma(self, ddataset2, data_path, rewrite_data_files):
        # make sure serialize looks correct
        session1 = pybmds.Session(dataset=ddataset2)
        session1.add_default_bayesian_models()
        session1.add_model_averaging()
        session1.execute_and_recommend()
        d = session1.to_dict()

        if rewrite_data_files:
            (data_path / "reports/session-dichotomous.json").write_text(
                session1.serialize().model_dump_json()
            )

        # spot check a few keys
        assert d["dataset"]["doses"] == [0.0, 50.0, 100.0, 150.0, 200.0]
        assert len(d["models"]) == 9
        assert list(d["models"][0].keys()) == ["name", "model_class", "settings", "results"]
        assert d["model_average"]["model_indexes"] == [0, 1, 2, 3, 4, 5, 6, 7, 8]
        assert "bmd" in d["model_average"]["results"]

        # ensure we can convert back to a session
        session2 = pybmds.Session.from_serialized(json.loads(json.dumps(d)))
        assert isinstance(session2, pybmds.Session)
        assert session2.dataset.doses == [0.0, 50.0, 100.0, 150.0, 200.0]
        assert len(session2.models) == 9
        assert session2.models[0].has_results is True

        # make sure we get the same result back after deserializing
        d1 = session1.serialize().model_dump()
        d2 = session2.serialize().model_dump()
        assert d1 == d2

        # df/docx
        df = session1.to_df()
        docx = session1.to_docx(session_inputs_table=True, all_models=True, bmd_cdf_table=True)

        if rewrite_data_files:
            df.to_excel(data_path / "reports/session-dichotomous-ma.xlsx", index=False)
            docx.save(data_path / "reports/session-dichotomous-ma.docx")

    # def test_dichotomous_ma_loud(self, ddataset2):
    #     session = pybmds.Session(dataset=ddataset2)
    #     session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)
    #     session.execute_and_recommend()
    #     assert session.model_average is not None

    def test_dichotomous_ma_rejects_mixed_priors(self, ddataset2):
        session = pybmds.Session(dataset=ddataset2)
        session.add_model(Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model(Models.Weibull, {"priors": PriorClass.bayesian_loud})

        with pytest.raises(ValueError, match="same prior_class|requires all models"):
            session.add_model_averaging()

    def test_continuous_ma(self, cdataset3, data_path, rewrite_data_files):
        # make sure serialize looks correct
        session1 = pybmds.Session(dataset=cdataset3)
        session1.add_default_bayesian_models()
        session1.add_model_averaging()
        session1.execute_and_recommend()
        d = session1.to_dict()

        if rewrite_data_files:
            (data_path / "reports/session-continuous.json").write_text(
                session1.serialize().model_dump_json()
            )

        # spot check a few keys
        assert d["dataset"]["doses"] == [0, 0.125, 0.25, 0.5, 1.0]
        assert len(d["models"]) == 10
        assert list(d["models"][0].keys()) == ["name", "model_class", "settings", "results"]
        assert d["model_average"]["model_indexes"] == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert "bmd" in d["model_average"]["results"]

        # ensure we can convert back to a session
        session2 = pybmds.Session.from_serialized(json.loads(json.dumps(d)))
        assert isinstance(session2, pybmds.Session)
        assert session2.dataset.doses == [0, 0.125, 0.25, 0.5, 1.0]
        assert len(session2.models) == 10
        assert session2.models[0].has_results is True

        # make sure we get the same result back after deserializing
        d1 = session1.serialize().model_dump()
        d2 = session2.serialize().model_dump()
        assert d1 == d2

        # df/docx
        df = session1.to_df()
        docx = session1.to_docx(session_inputs_table=True, all_models=True, bmd_cdf_table=False)

        if rewrite_data_files:
            df.to_excel(data_path / "reports/session-continuous-ma.xlsx", index=False)
            docx.save(data_path / "reports/session-continuous-ma.docx")

    def test_continuous_ma_rejects_bayesian_priors(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)

        # Add eligible continuous models with "regular" bayesian priors
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian}
        )
        session.add_model(Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian})

        with pytest.raises(ValueError, match=r"Continuous model averaging requires.*bayesian_loud"):
            session.add_model_averaging()

    def test_continuous_ma_allows_bayesian_loud_priors(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)

        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )

        session.add_model_averaging()
        session.execute_and_recommend()
        d = session.to_dict()

        assert session.model_average is not None
        assert "bmd" in d["model_average"]["results"]

    def test_continuous_ma_manual_models_keep_hill_with_efsa_models(self, cdataset3, monkeypatch):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Power, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.LMS2, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.InverseExponential,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )
        session.add_model_averaging()
        session.execute()

        assert [model.name() for model in session.model_average.models] == [
            "Power (CV)",
            "Hill (CV)",
            "Power (NCV)",
            "Hill (NCV)",
            "LMS 2-Stage (NCV)",
            "Inverse Exponential (NCV)",
        ]
        assert len(session.model_average.results.priors) == len(session.model_average.models)
        assert len(session.model_average.results.posteriors) == len(session.model_average.models)
        idata = pybmds.plotting.LOUD.model_average_to_inferencedata(session)
        assert list(idata.posterior.coords["model"].values) == [
            "Power (CV)",
            "Hill (CV)",
            "Power (NCV)",
            "Hill (NCV)",
            "LMS 2-Stage (NCV)",
            "Inverse Exponential (NCV)",
        ]

        def fake_get_model_average_figures(_session, n_chains=1):
            assert _session is session
            assert n_chains == 1

            def fig():
                return plt.figure()

            return {
                "posterior": fig(),
                "overlay": fig(),
                "bmd_summary": pd.DataFrame({"median": [1.23]}, index=["MA_BMD"]),
                "parameter_groups": [],
                "alpha": 0.05,
                "hdi_prob": 0.9,
            }

        monkeypatch.setattr(
            pybmds.session, "get_model_average_figures", fake_get_model_average_figures
        )

        docx = session.to_docx(citation=False)
        bayesian_table = docx.tables[1]

        assert bayesian_table.cell(1, 1).text != "-"
        assert bayesian_table.cell(2, 1).text != "-"
        assert bayesian_table.cell(3, 1).text != "-"
        assert bayesian_table.cell(4, 1).text != "-"
        assert bayesian_table.cell(5, 1).text != "-"
        assert bayesian_table.cell(6, 1).text != "-"

    def test_continuous_ma_syncs_loud_per_model_results_back_to_models(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.LMS2, {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.InverseExponential,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )
        session.add_model_averaging()
        session.execute()

        loud_lookup = {
            model.name(): session.model_average.results.model_summary(idx, model.settings.alpha)
            for idx, model in enumerate(session.model_average.models)
        }

        for model in session.models:
            assert model.has_results is True
            assert "Model has not successfully executed" not in model.text()
            assert "Model Parameters:" in model.text()

            summary = loud_lookup[model.name()]
            assert model.results.bmdl == pytest.approx(summary.bmdl)
            assert model.results.bmd == pytest.approx(summary.bmd)
            assert model.results.bmdu == pytest.approx(summary.bmdu)
            assert len(model.results.fit.bmd_dist) > 0
            assert len(model.results.parameters.names) == len(model.results.parameters.values)
            assert len(model.results.parameters.names) == len(model.results.parameters.se)
            assert np.isfinite(model.results.plotting.dr_y).all()
            assert np.unique(model.results.plotting.dr_y).size > 1

    def test_continuous_ma_default_include_efsa_excludes_existing_hill(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )

        session.add_default_bayesian_models(
            include_efsa=True,
            model_average=False,
            prior_class=PriorClass.bayesian_loud,
        )
        session.add_model_averaging()

        assert "Hill (CV)" in [model.name() for model in session.models]
        assert "Hill (CV)" not in [model.name() for model in session.model_average.models]

    def test_to_docx_all_models_skips_cdf_for_unsuccessful_models(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_default_bayesian_models(
            include_efsa=True,
            model_average=False,
            prior_class=PriorClass.bayesian_loud,
        )
        session.execute()

        assert any(not model.has_results for model in session.models)

        docx = session.to_docx(citation=False, all_models=True, bmd_cdf_table=True)

        assert docx is not None

    def test_continuous_loud_docx_reporting_spacing(self, cdataset3, monkeypatch):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        def fake_get_model_average_figures(_session, n_chains=1):
            assert _session is session
            assert n_chains == 1

            def fig():
                return plt.figure()

            return {
                "posterior": fig(),
                "overlay": fig(),
                "bmd_summary": pd.DataFrame(
                    {
                        "median": [1.23],
                        "eti_5%": [1.0],
                        "eti_95%": [1.5],
                        "r_hat": [1.0],
                        "ess_bulk": [90.0],
                        "ess_tail": [80.0],
                    },
                    index=["MA_BMD"],
                ),
                "parameter_groups": [
                    {
                        "name": "Power",
                        "summary": pd.DataFrame(
                            {
                                "Model": ["CV", "NCV"],
                                "Parameter": ["g", "rho"],
                                "median": [1.23, 2.34],
                                "r_hat": [1.0, 1.0],
                                "ess_bulk": [90.0, 88.0],
                                "ess_tail": [80.0, 79.0],
                            }
                        ),
                        "trace_figure": fig(),
                    }
                ],
                "alpha": 0.05,
                "hdi_prob": 0.9,
            }

        monkeypatch.setattr(
            pybmds.session, "get_model_average_figures", fake_get_model_average_figures
        )

        docx = session.to_docx(citation=False, bmd_cdf_table=False)
        spacing_by_text = {
            paragraph.text: paragraph.paragraph_format.space_before.pt
            if paragraph.paragraph_format.space_before is not None
            else None
            for paragraph in docx.paragraphs
        }

        assert spacing_by_text["Model Averaging Diagnostics (LOUD)"] == 6.0
        assert spacing_by_text["Posterior distribution of model-averaged BMD"] == 6.0
        assert (
            spacing_by_text["Overlay of model-specific and model-averaged BMD distributions"] == 6.0
        )
        assert spacing_by_text["Summary statistics for BMD and model-averaged BMD"] == 6.0
        assert spacing_by_text["Power model parameters"] == 6.0
        assert "Power model parameter visualizations" not in spacing_by_text

    def test_continuous_ma_to_df_handles_models_without_toi_tests(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.ExponentialM3,
            {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.MultiplicativeHill,
            {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.ExponentialM3,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.MultiplicativeHill,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.LMS2,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )
        session.add_model_averaging()
        session.execute()

        df = session.to_df()

        blank = pybmds.constants.BMDS_BLANK_VALUE
        mh_row = df[df.model_name.str.startswith("Multiplicative Hill")].iloc[0]
        lms_row = df[df.model_name.str.startswith("LMS 2-Stage")].iloc[0]

        assert mh_row.p_value1 == blank
        assert mh_row.p_value4 == blank
        assert mh_row.model_dof == blank
        assert lms_row.p_value1 == blank
        assert lms_row.p_value4 == blank
        assert lms_row.model_dof == blank

    def test_to_docx_can_include_grouped_parameter_visualizations(self, monkeypatch, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        def fake_get_model_average_figures(_session, n_chains=1):
            def fig():
                return plt.figure()

            return {
                "posterior": fig(),
                "overlay": fig(),
                "bmd_summary": pd.DataFrame({"median": [1.23]}, index=["MA_BMD"]),
                "parameter_groups": [
                    {
                        "name": "Power",
                        "summary": pd.DataFrame(
                            {
                                "Model": ["CV"],
                                "Parameter": ["g"],
                                "median": [1.23],
                                "r_hat": [1.0],
                                "ess_bulk": [90.0],
                                "ess_tail": [80.0],
                            }
                        ),
                        "trace_figure": fig(),
                    }
                ],
                "alpha": 0.05,
                "hdi_prob": 0.9,
            }

        monkeypatch.setattr(
            pybmds.session, "get_model_average_figures", fake_get_model_average_figures
        )

        docx = session.to_docx(citation=False, parameter_visualizations=True)
        paragraph_text = [paragraph.text for paragraph in docx.paragraphs]

        assert "Power model parameters" in paragraph_text
        assert "Power model parameter visualizations" in paragraph_text

    def test_to_docx_can_hide_parameter_tables(self, monkeypatch, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()

        def fake_get_model_average_figures(_session, n_chains=1):
            def fig():
                return plt.figure()

            return {
                "posterior": fig(),
                "overlay": fig(),
                "bmd_summary": pd.DataFrame({"median": [1.23], "r_hat": [1.0]}, index=["MA_BMD"]),
                "parameter_groups": [
                    {
                        "name": "Power",
                        "summary": pd.DataFrame({"Model": ["CV"], "Parameter": ["g"]}),
                        "trace_figure": fig(),
                    }
                ],
                "alpha": 0.05,
                "hdi_prob": 0.9,
            }

        monkeypatch.setattr(
            pybmds.session, "get_model_average_figures", fake_get_model_average_figures
        )

        docx = session.to_docx(citation=False, parameter_tables=False)
        paragraph_text = [paragraph.text for paragraph in docx.paragraphs]

        assert "Power model parameters" not in paragraph_text

    def test_continuous_manual_models_use_disttype_in_default_names(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.Power,
            {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.Hill,
            {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.Power,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )
        session.add_model(
            pybmds.Models.Hill,
            {"disttype": DistType.normal_ncv, "priors": PriorClass.bayesian_loud},
        )

        assert [model.name() for model in session.models] == [
            "Power (CV)",
            "Hill (CV)",
            "Power (NCV)",
            "Hill (NCV)",
        ]

    def test_continuous_cma_with_efsa(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)

        session.add_default_bayesian_models(include_efsa=True)
        assert getattr(session, "models", None) is not None
        assert len(session.models) > 0, "No models were added to the session"

        session.add_model_averaging()

        ma = session.model_average
        assert session.model_average is not None, "model_average was not created"
        ma = session.model_average

        assert getattr(ma, "models", None) is not None, "model_average.models is missing"
        assert len(ma.models) > 0, "model_average has no models"

        from pybmds.constants import DistType

        seen = {}
        for m in ma.models:
            seen.setdefault(m.__class__.__name__, set()).add(m.settings.disttype)

        # Spot-check one EFSA model
        assert "MultiplicativeHill" in seen, f"MultiplicativeHill not found; saw: {sorted(seen)}"
        assert seen["MultiplicativeHill"] == {
            DistType.normal,
            DistType.normal_ncv,
            DistType.log_normal,
        }, f"MultiplicativeHill disttypes wrong: {seen['MultiplicativeHill']}"

        # BMDS models must still be present
        assert "Power" in seen, f"Power not found; saw: {sorted(seen)}"
        assert "ExponentialM3" in seen, f"ExponentialM3 not found; saw: {sorted(seen)}"

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
        assert int(version.split(".")[0]) >= 24  # assume dll in format "YY.MM..."


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
    def test_continuous_ma_colorize(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            pybmds.Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        session.execute()
        return session.plot(colorize=True)

    @pytest.mark.mpl_image_compare
    def test_continuous_ma(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            pybmds.Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
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
