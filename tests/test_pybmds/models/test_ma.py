import json
from types import SimpleNamespace

import numpy as np
import pytest

import pybmds
import pybmds.models.ma as ma_module
from pybmds.constants import PriorClass
from pybmds.models.ma import BmdModelAveragingDichotomous
from pybmds.types.dichotomous import DichotomousModelSettings
from pybmds.types.ma import DichotomousModelAverage, DichotomousModelAverageResult


def loud_dma_dataset():
    return pybmds.DichotomousDataset(
        doses=[0, 0.25, 0.75, 0.85, 1],
        ns=[20, 20, 20, 20, 20],
        incidences=[0, 1, 7, 15, 19],
    )


class TestDichotomousMa:
    def test_model_settings_accept_defaults_instances_and_dicts(self):
        defaults = BmdModelAveragingDichotomous.get_model_settings(None, None)
        assert isinstance(defaults, DichotomousModelSettings)
        assert BmdModelAveragingDichotomous.get_model_settings(None, defaults) is defaults
        assert BmdModelAveragingDichotomous.get_model_settings(None, {"alpha": 0.1}).alpha == 0.1

    def test_default_dichotomous_loud_session_and_ma_models_use_same_order(self, ddataset2):
        session = pybmds.Session(dataset=ddataset2)
        session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)

        expected = [
            "Hill",
            "Gamma",
            "Logistic",
            "LogLogistic",
            "LogProbit",
            "Multistage 2",
            "Probit",
            "Quantal Linear",
            "Weibull",
        ]
        assert [model.name() for model in session.models] == expected
        assert [model.name() for model in session.model_average.models] == expected

    def test_dichotomous_loud_ma_uses_fixed_model_order(self, ddataset2):
        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Weibull, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.DichotomousHill, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian_loud})
        session.add_model_averaging(weights=[0.2, 0.3, 0.1, 0.4])

        assert [model.name() for model in session.models] == [
            "Weibull",
            "Logistic",
            "Hill",
            "Probit",
        ]
        assert [model.name() for model in session.model_average.models] == [
            "Hill",
            "Logistic",
            "Probit",
            "Weibull",
        ]
        assert np.allclose(session.ma_weights, [0.2, 0.3, 0.1, 0.4])

    def test_dichotomous_loud_ma_remaps_weights_to_fixed_model_order(self, ddataset2, monkeypatch):
        captured = {}

        class FakeDichotomousModelAverage:
            def __init__(self, dataset, models, model_weights, prior_class, weight_option):
                captured["model_names"] = [model.name() for model in models]
                captured["model_weights"] = np.asarray(model_weights, dtype=float)
                self.result = SimpleNamespace(models=[SimpleNamespace() for _ in models])

            def execute(self):
                return self

        fake_result = SimpleNamespace(sync_model_result=lambda *args: None)
        monkeypatch.setattr(ma_module, "DichotomousModelAverage", FakeDichotomousModelAverage)
        monkeypatch.setattr(
            ma_module.DichotomousModelAverageResult,
            "from_cpp",
            classmethod(lambda cls, structs, model_results: fake_result),
        )

        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Weibull, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.DichotomousHill, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian_loud})
        session.add_model_averaging(weights=[0.2, 0.3, 0.1, 0.4])

        session.model_average.execute()

        assert captured["model_names"] == ["Hill", "Logistic", "Probit", "Weibull"]
        assert np.allclose(captured["model_weights"], [0.1, 0.3, 0.4, 0.2])

    def test_regular_dichotomous_bayesian_ma_keeps_session_order_and_model_priors(
        self, ddataset2, monkeypatch
    ):
        captured = {}

        class FakeDichotomousModelAverage:
            def __init__(self, dataset, models, model_weights, prior_class, weight_option):
                captured["model_names"] = [model.name() for model in models]
                captured["model_weights"] = np.asarray(model_weights, dtype=float)
                self.result = SimpleNamespace(models=[SimpleNamespace() for _ in models])

            def execute(self):
                return self

        monkeypatch.setattr(ma_module, "DichotomousModelAverage", FakeDichotomousModelAverage)
        monkeypatch.setattr(
            ma_module.DichotomousModelAverageResult,
            "from_cpp",
            classmethod(lambda cls, structs, model_results: SimpleNamespace()),
        )

        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Weibull, {"priors": PriorClass.bayesian})
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian})
        session.add_model_averaging(weights=[0.2, 0.3, 0.5])

        session.model_average.execute()

        assert [model.name() for model in session.model_average.models] == [
            "Weibull",
            "Logistic",
            "Probit",
        ]
        assert captured["model_names"] == ["Weibull", "Logistic", "Probit"]
        assert np.allclose(captured["model_weights"], [0.2, 0.3, 0.5])

    def test_dichotomous_ma_session(self, ddataset2):
        # check execution and it can be json serialized
        session = pybmds.Session(dataset=ddataset2)
        session.add_default_bayesian_models()
        session.execute()
        d = session.to_dict()
        assert isinstance(json.dumps(d), str)

        # check bmd values exist and are valid
        res = session.model_average.results
        assert np.allclose([57.1, 65.9, 75.0], [res.bmdl, res.bmd, res.bmdu], atol=5)

    def test_prior_weights(self, ddataset2):
        # default; equal weights
        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian})
        session.add_model_averaging()
        assert np.allclose(session.ma_weights, [0.5, 0.5])
        session.execute()
        assert np.allclose(session.model_average.results.priors, [0.5, 0.5])
        assert np.allclose(session.model_average.results.posteriors, [0.11, 0.89], atol=0.05)

        # custom; propagate through results
        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian})
        session.add_model_averaging(weights=[0.9, 0.1])
        assert np.allclose(session.ma_weights, [0.9, 0.1])
        session.execute()
        assert np.allclose(session.model_average.results.priors, [0.9, 0.1])
        assert np.allclose(session.model_average.results.posteriors, [0.53, 0.47], atol=0.05)

    def test_dichotomous_loud_bridge_weights_are_reported(self):
        dataset = pybmds.DichotomousDataset(
            doses=[0, 0.25, 0.5, 0.75, 1],
            ns=[20, 20, 20, 20, 20],
            incidences=[0, 3, 10, 15, 19],
        )
        session = pybmds.Session(dataset=dataset)
        session.add_default_bayesian_models(
            prior_class=PriorClass.bayesian_loud, weight_option="bridge"
        )
        session.execute()

        posteriors = session.model_average.results.posteriors
        assert np.isfinite(posteriors).all()
        assert posteriors.sum() == pytest.approx(1)

        docx = session.to_docx(citation=False)
        bayesian_table = next(
            table
            for table in docx.tables
            if table.rows[0].cells[0].text == "Model"
            and table.rows[0].cells[2].text == "Posterior Weights"
        )
        for row in bayesian_table.rows[1:-1]:
            assert row.cells[2].text != "-"

    def test_dichotomous_ma_datatype_matches_prior_class(self, ddataset2):
        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian})
        session.add_model_averaging()
        for model in session.models:
            model.execute_job()
        session.model_average.structs = DichotomousModelAverage(
            session.dataset,
            session.model_average.models,
            session.ma_weights,
            session.model_average.settings.priors.prior_class,
        )
        assert session.model_average.structs.average.datatype == 0

        session = pybmds.Session(dataset=ddataset2)
        session.add_model(pybmds.Models.Logistic, {"priors": PriorClass.bayesian_loud})
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian_loud})
        session.add_model_averaging()
        for model in session.models:
            model.execute_job()
        session.model_average.structs = DichotomousModelAverage(
            session.dataset,
            session.model_average.models,
            session.ma_weights,
            session.model_average.settings.priors.prior_class,
        )
        assert session.model_average.structs.average.datatype == 4

    def test_loud_mcmc_settings_are_staged_for_dichotomous_ma(self, ddataset2):
        session = pybmds.Session(dataset=ddataset2)
        session.add_model(
            pybmds.Models.Logistic,
            {"priors": PriorClass.bayesian_loud, "n_chains": 4, "samples": 25_000, "seed": 123},
        )
        session.add_model(pybmds.Models.Probit, {"priors": PriorClass.bayesian_loud})
        session.add_model_averaging()
        for model in session.models:
            model.execute_job()

        structs = DichotomousModelAverage(
            session.dataset,
            session.model_average.models,
            session.ma_weights,
            session.model_average.settings.priors.prior_class,
        )

        assert structs.n_chains == 4
        assert structs.seed == 123
        assert structs.analysis.samples == 25_000
        assert structs.analysis.chains == 4

    def test_dichotomous_loud_ma_session_syncs_results(self):
        session = pybmds.Session(dataset=loud_dma_dataset())
        session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)

        session.execute()

        ma = session.model_average
        assert ma is not None
        assert ma.results is not None
        assert ma.results.bmd > 0
        assert ma.results.bmdl > 0
        assert len(ma.results.bmd_dist) > 0
        assert len(ma.results.model_bmd_dist) == len(ma.models)
        assert len(ma.results.model_parm_dist) == len(ma.models)

        summary = ma.results.model_summary(0, ma.models[0].settings.alpha)
        assert summary.bmdl > 0
        assert summary.bmd > 0
        assert summary.prior == pytest.approx(ma.results.priors[0])
        assert summary.posterior == pytest.approx(ma.results.posteriors[0])

        summary_lookup = {
            model: ma.results.model_summary(idx, model.settings.alpha)
            for idx, model in enumerate(ma.models)
        }
        for model in session.models:
            assert model.has_results is True
            summary = summary_lookup[model]
            assert model.results.bmd == pytest.approx(summary.bmd)
            if summary.bmd > 0:
                assert len(model.results.fit.bmd_dist) > 0
            assert len(model.results.parameters.names) == len(model.results.parameters.values)
            assert len(model.results.parameters.names) == len(model.results.parameters.se)

        for idx, model in enumerate(ma.models):
            assert ma.results.model_bmd_dist[idx].shape[0] == len(ma.results.bmd_dist)
            assert ma.results.model_parm_dist[idx].shape[0] == len(ma.results.bmd_dist)
            assert ma.results.model_parm_dist[idx].shape[1] == len(model.results.parameters.names)

    def test_dichotomous_loud_ma_docx_uses_real_bayesian_table(self):
        session = pybmds.Session(dataset=loud_dma_dataset())
        session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)
        session.execute()

        docx = session.to_docx(citation=False, bmd_cdf_table=False)
        paragraph_text = [paragraph.text for paragraph in docx.paragraphs]
        parameter_model_labels = [
            row.cells[0].text
            for table in docx.tables
            if len(table.rows) > 1
            and [cell.text for cell in table.rows[0].cells[:2]] == ["Model", "Parameter"]
            for row in table.rows[1:]
        ]

        assert "Model Averaging Diagnostics (LOUD)" in paragraph_text
        assert len(docx.tables) > 0
        assert "None" not in parameter_model_labels
        assert any(label in parameter_model_labels for label in ["Logistic", "Weibull"])

    def test_dichotomous_loud_ma_all_models_report_uses_loud_model_results(self):
        session = pybmds.Session(dataset=loud_dma_dataset())
        session.add_default_bayesian_models(prior_class=PriorClass.bayesian_loud)
        session.execute()

        docx = session.to_docx(citation=False, all_models=True, bmd_cdf_table=True)
        paragraph_text = "\n".join(paragraph.text for paragraph in docx.paragraphs)

        assert "Individual Model Results" in paragraph_text
        assert "LOUD Model-Average Weights" in paragraph_text

        ma = session.model_average
        summary_lookup = {
            model: ma.results.model_summary(idx, model.settings.alpha)
            for idx, model in enumerate(ma.models)
        }
        p_value_lookup = {
            model: ma.results.model_p_values[idx] for idx, model in enumerate(ma.models)
        }

        for model in session.models:
            summary = summary_lookup[model]
            assert model.results.bmdl == pytest.approx(summary.bmdl)
            assert model.results.bmd == pytest.approx(summary.bmd)
            assert model.results.bmdu == pytest.approx(summary.bmdu)
            assert model.results.gof.p_value == pytest.approx(p_value_lookup[model])
            assert len(model.results.fit.bmd_dist) > 0

    def test_loud_ma_serializes_finite_draws_for_json(self):
        result = DichotomousModelAverageResult(
            bmdl=0.1,
            bmd=0.2,
            bmdu=0.3,
            bmdl_y=0.01,
            bmd_y=0.02,
            bmdu_y=0.03,
            bmd_dist=np.array([0.1, np.inf, 0.3]),
            priors=np.array([1.0]),
            posteriors=np.array([1.0]),
            model_bmd_dist=[np.array([0.1, np.inf, 0.3])],
            model_parm_dist=[
                np.array(
                    [
                        [1.0, 2.0],
                        [np.nan, 3.0],
                        [4.0, np.inf],
                    ]
                )
            ],
            model_p_values=[0.5],
            dr_x=np.array([0.0, 1.0]),
            dr_y=np.array([0.0, 1.0]),
        )

        assert np.isinf(result.bmd_dist[1])
        assert np.isinf(result.model_bmd_dist[0][1])
        assert np.isnan(result.model_parm_dist[0][1, 0])

        payload = result.model_dump()
        assert payload["bmd_dist"] == [0.1, None, 0.3]
        assert payload["model_bmd_dist"] == [[0.1, None, 0.3]]
        assert payload["model_parm_dist"] == [[[1.0, 2.0], [None, 3.0], [4.0, None]]]
        assert isinstance(json.dumps(payload, allow_nan=False), str)

        rehydrated = DichotomousModelAverageResult.model_validate(payload)
        assert len(rehydrated.bmd_dist) == 3
        assert len(rehydrated.model_bmd_dist[0]) == 3
        assert rehydrated.model_parm_dist[0].shape == (3, 2)
        assert np.isnan(rehydrated.bmd_dist[1])
        assert np.isnan(rehydrated.model_bmd_dist[0][1])
        assert np.isnan(rehydrated.model_parm_dist[0][1, 0])

    def test_loud_helper_branches(self):
        assert DichotomousModelAverageResult._param_names(
            SimpleNamespace(get_param_names=lambda: ["a"]), 3
        ) == ["a", "param_2", "param_3"]

        class BadCurve:
            def dr_curve(self, *_args):
                raise ValueError("bad")

        assert DichotomousModelAverageResult._safe_dr_curve(BadCurve(), np.array([1.0]), []) is None

        one_d = DichotomousModelAverageResult._resize_draws(np.array([1.0]), 3)
        np.testing.assert_array_equal(one_d[:1], [1.0])
        assert np.isnan(one_d[1])
        two_d = DichotomousModelAverageResult._resize_draws(np.array([[1.0, 2.0]]), 2)
        assert two_d.shape == (2, 2)

        bmd, parms = DichotomousModelAverageResult._apply_paired_valid_draw_mask(
            np.array([0.1, 0.2, 0.3]),
            np.array([[1.0, 2.0], [np.nan, 3.0], [4.0, np.inf]]),
        )
        np.testing.assert_array_equal(bmd, np.array([0.1, np.nan, np.nan]))
        np.testing.assert_array_equal(
            parms,
            np.array([[1.0, 2.0], [np.nan, np.nan], [np.nan, np.nan]]),
        )

        np.testing.assert_array_equal(
            DichotomousModelAverageResult._loud_weights([np.nan, -9999.0]),
            np.array([-9999.0, -9999.0]),
        )
        np.testing.assert_array_equal(
            DichotomousModelAverageResult._loud_weights([1000.0, np.inf]),
            np.array([1.0, 0.0]),
        )

        existing = np.array([0.25, 0.75])
        np.testing.assert_array_equal(
            DichotomousModelAverageResult._loud_posteriors(existing, 1, [1.0, 2.0], [3.0, 4.0]),
            existing,
        )
        assert (
            DichotomousModelAverageResult._loud_posteriors(
                np.array([np.nan, np.nan]), 2, [1.0, 2.0], [3.0, 4.0]
            )[1]
            > 0
        )
        assert np.all(
            DichotomousModelAverageResult._loud_posteriors(
                np.array([np.nan, np.nan]), 3, [np.nan, -9999.0], [np.nan, -9999.0]
            )
            == -9999.0
        )

    def test_loud_draw_helpers_cover_chain_and_fallback_shapes(self):
        chains = [
            SimpleNamespace(BMD=np.array([1.0, 2.0]), parms=np.array([1.0, 2.0, 3.0, 4.0])),
            SimpleNamespace(BMD=np.array([3.0, 4.0]), parms=np.array([[5.0, 6.0], [7.0, 8.0]])),
        ]

        bmd, parms = DichotomousModelAverageResult._loud_draws(chains, n_chains=2)

        assert bmd.shape == (2, 2)
        assert parms.shape == (2, 2, 2)

        loud = SimpleNamespace(
            BMD=np.array([1.0, 2.0, 3.0, 4.0]),
            parms=np.array([[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]),
        )
        bmd, parms = DichotomousModelAverageResult._loud_draws(loud, n_chains=2)
        assert bmd.shape == (2, 2)
        assert parms.shape == (4, 2)

        combined = SimpleNamespace(BMD=[1.0])
        result = SimpleNamespace(combinedLoudRes=combined, loudRes=chains)
        assert DichotomousModelAverageResult._combined_loud_result(result) is combined
