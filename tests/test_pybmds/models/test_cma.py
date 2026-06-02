import json
from types import SimpleNamespace

import numpy as np

import pybmds
from pybmds.constants import DistType, PriorClass
from pybmds.types.cma import ContinuousModelAverage, ContinuousModelAverageResult


## TO DO - to change when we have actual results from models
class TestContinuousMa:
    def test_continuous_individual_loud_ma_uses_individual_data(self, cidataset):
        session = pybmds.Session(dataset=cidataset)
        session.add_model(pybmds.Models.Power, {"priors": PriorClass.bayesian_loud})
        session.models[0].settings.samples = 20
        session.models[0].settings.burnin = 5
        session.models[0].prepare_for_loud_model_average()
        session.add_model_averaging()

        structs = ContinuousModelAverage(
            session.dataset,
            session.model_average.models,
            session.ma_weights,
            session.weight_option,
        )

        assert structs.analysis.n == len(cidataset.individual_doses)
        assert structs.analysis.doses == cidataset.individual_doses
        assert structs.analysis.Y == cidataset.responses
        assert structs.analysis.sd == []
        assert structs.analysis.n_group == []

    def test_loud_mcmc_settings_are_staged_for_continuous_ma(self, cdataset3):
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.Power,
            {
                "disttype": DistType.normal,
                "priors": PriorClass.bayesian_loud,
                "n_chains": 4,
                "seed": 123,
            },
        )
        session.models[0].prepare_for_loud_model_average()
        session.add_model_averaging()

        structs = ContinuousModelAverage(
            session.dataset,
            session.model_average.models,
            session.ma_weights,
            session.weight_option,
        )

        assert structs.n_chains == 4
        assert structs.seed == 123
        assert structs.average.seed == 123

    def test_continuous_ma_session(self, cdataset3):
        # check execution and it can be json serialized
        session = pybmds.Session(dataset=cdataset3)
        session.add_default_bayesian_models()
        session.execute()
        d = session.to_dict()
        assert isinstance(json.dumps(d), str)

        # check bmd values exist and are valid
        res = session.model_average.results
        assert np.allclose([0.036, 0.079, 0.22], [res.bmdl, res.bmd, res.bmdu], atol=5)

    def test_prior_weights(self, cdataset3):
        # default; equal weights
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            pybmds.Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging()
        assert np.allclose(session.ma_weights, [0.5, 0.5])
        session.execute()
        assert np.allclose(session.model_average.results.priors, [0.5, 0.5])
        assert np.allclose(session.model_average.results.posteriors, [0.0043, 0.9957], atol=0.05)

        # custom; propagate through results
        session = pybmds.Session(dataset=cdataset3)
        session.add_model(
            pybmds.Models.Power, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model(
            pybmds.Models.Hill, {"disttype": DistType.normal, "priors": PriorClass.bayesian_loud}
        )
        session.add_model_averaging(weights=[0.9, 0.1])
        assert np.allclose(session.ma_weights, [0.9, 0.1])
        session.execute()
        assert np.allclose(session.model_average.results.priors, [0.9, 0.1])
        assert np.allclose(session.model_average.results.posteriors, [0.0043, 0.9957], atol=0.05)

    def test_loud_ma_preserves_nonfinite_parameter_draw_rows(self):
        analysis = SimpleNamespace(
            result=SimpleNamespace(
                bmd_dist=np.array([0.1, np.inf, 0.3]),
                post_probs=np.array([1.0]),
                bmdsRes=SimpleNamespace(BMDL_MA=0.1, BMD_MA=0.2, BMDU_MA=0.3),
                models=[
                    SimpleNamespace(
                        loudRes=SimpleNamespace(
                            BMD=np.array([0.1, np.inf, 0.3]),
                            parms=np.array(
                                [
                                    [1.0, 2.0],
                                    [np.nan, 3.0],
                                    [4.0, np.inf],
                                ]
                            ),
                            pval=0.5,
                        )
                    )
                ],
            ),
            average=SimpleNamespace(modelPriors=np.array([1.0])),
        )
        model_results = [
            SimpleNamespace(
                plotting=SimpleNamespace(
                    dr_x=np.array([0.0, 1.0]),
                    dr_y=np.array([10.0, 9.0]),
                )
            )
        ]

        result = ContinuousModelAverageResult.from_cpp(analysis, model_results)

        assert result.model_parm_dist[0].shape == (3, 2)
        assert np.isnan(result.model_parm_dist[0][1, 0])
        assert np.isinf(result.model_parm_dist[0][2, 1])
        assert np.isinf(result.bmd_dist[1])
        assert np.isinf(result.model_bmd_dist[0][1])

        payload = result.model_dump()
        assert payload["bmd_dist"] == [0.1, None, 0.3]
        assert payload["model_bmd_dist"] == [[0.1, None, 0.3]]
        assert payload["model_parm_dist"] == [[[1.0, 2.0], [None, 3.0], [4.0, None]]]
        assert isinstance(json.dumps(payload, allow_nan=False), str)

        rehydrated = ContinuousModelAverageResult.model_validate(payload)
        assert len(rehydrated.bmd_dist) == 3
        assert len(rehydrated.model_bmd_dist[0]) == 3
        assert rehydrated.model_parm_dist[0].shape == (3, 2)
        assert np.isnan(rehydrated.bmd_dist[1])
        assert np.isnan(rehydrated.model_bmd_dist[0][1])
        assert np.isnan(rehydrated.model_parm_dist[0][1, 0])
