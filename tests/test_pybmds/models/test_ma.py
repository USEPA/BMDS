import json

import numpy as np
import pytest

import pybmds
from pybmds.constants import PriorClass
from pybmds.types.ma import DichotomousModelAverage


def loud_dma_dataset():
    return pybmds.DichotomousDataset(
        doses=[0, 0.25, 0.75, 0.85, 1],
        ns=[20, 20, 20, 20, 20],
        incidences=[0, 1, 7, 15, 19],
    )


class TestDichotomousMa:
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

        summary = ma.results.model_summary(0, session.models[0].settings.alpha)
        assert summary.bmdl > 0
        assert summary.bmd > 0
        assert summary.prior == pytest.approx(ma.results.priors[0])
        assert summary.posterior == pytest.approx(ma.results.posteriors[0])

        for idx, model in enumerate(session.models):
            assert model.has_results is True
            summary = ma.results.model_summary(idx, model.settings.alpha)
            assert model.results.bmd == pytest.approx(summary.bmd)
            if summary.bmd > 0:
                assert len(model.results.fit.bmd_dist) > 0
            assert len(model.results.parameters.names) == len(model.results.parameters.values)
            assert len(model.results.parameters.names) == len(model.results.parameters.se)
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
