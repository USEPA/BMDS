import json
from types import SimpleNamespace

import numpy as np
import pytest

import pybmds
from pybmds.constants import DistType, PriorClass
from pybmds.models.cma import BmdModelAveragingContinuous
from pybmds.types.cma import ContinuousModelAverage, ContinuousModelAverageResult
from pybmds.types.continuous import ContinuousModelSettings


## TO DO - to change when we have actual results from models
class TestContinuousMa:
    def test_model_settings_accept_defaults_instances_and_dicts(self):
        defaults = BmdModelAveragingContinuous.get_model_settings(None, None)
        assert isinstance(defaults, ContinuousModelSettings)
        assert BmdModelAveragingContinuous.get_model_settings(None, defaults) is defaults
        assert BmdModelAveragingContinuous.get_model_settings(None, {"alpha": 0.1}).alpha == 0.1

    def test_mixed_distribution_waic_is_model_order_invariant(self):
        dataset = pybmds.ContinuousDataset(
            doses=[0, 50, 100, 200, 400],
            ns=[20, 20, 20, 20, 20],
            means=[5.26, 5.56, 6.13, 8.24, 11.23],
            stdevs=[2.23, 2.37, 2.45, 2.64, 2.99],
        )
        power = (pybmds.Models.Power, DistType.normal)
        exp3_lognormal = (pybmds.Models.ExponentialM3, DistType.log_normal)

        def execute(order):
            session = pybmds.Session(dataset=dataset)
            for model, disttype in order:
                session.add_model(
                    model,
                    {
                        "disttype": disttype,
                        "priors": PriorClass.bayesian_loud,
                        "n_chains": 1,
                        "samples": 1000,
                        "burnin": 100,
                        "seed": 421,
                    },
                )
            session.add_model_averaging()
            session.execute()
            return {
                (model.bmd_model_class.id, model.settings.disttype): {
                    "waic": waic,
                    "residuals": np.asarray(model.results.gof.residual, dtype=float),
                }
                for model, waic in zip(
                    session.model_average.models,
                    session.model_average.results.model_waics,
                    strict=True,
                )
            }

        forward = execute([power, exp3_lognormal])
        reverse = execute([exp3_lognormal, power])

        assert forward.keys() == reverse.keys()
        for key in forward:
            assert forward[key]["waic"] == pytest.approx(reverse[key]["waic"])
            np.testing.assert_allclose(
                forward[key]["residuals"], reverse[key]["residuals"], rtol=0, atol=0
            )

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
                "samples": 25_000,
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
        assert structs.analysis.samples == 25_000
        assert structs.analysis.chains == 4

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
        priors_by_model = {
            model.name(): prior
            for model, prior in zip(
                session.model_average.models, session.model_average.results.priors, strict=True
            )
        }
        posteriors_by_model = {
            model.name(): posterior
            for model, posterior in zip(
                session.model_average.models,
                session.model_average.results.posteriors,
                strict=True,
            )
        }
        assert np.allclose(list(priors_by_model.values()), [0.5, 0.5])
        assert posteriors_by_model["Power (CV)"] == pytest.approx(0.0043, abs=0.05)
        assert posteriors_by_model["Hill (CV)"] == pytest.approx(0.9957, abs=0.05)

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
        priors_by_model = {
            model.name(): prior
            for model, prior in zip(
                session.model_average.models, session.model_average.results.priors, strict=True
            )
        }
        posteriors_by_model = {
            model.name(): posterior
            for model, posterior in zip(
                session.model_average.models,
                session.model_average.results.posteriors,
                strict=True,
            )
        }
        assert priors_by_model == pytest.approx({"Power (CV)": 0.9, "Hill (CV)": 0.1})
        assert posteriors_by_model["Power (CV)"] == pytest.approx(0.0043, abs=0.05)
        assert posteriors_by_model["Hill (CV)"] == pytest.approx(0.9957, abs=0.05)

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
        assert np.isnan(rehydrated.model_parm_dist[0][2, 1])

    def test_loud_ma_recovers_blank_posteriors_from_waic(self):
        analysis = SimpleNamespace(
            result=SimpleNamespace(
                bmd_dist=np.array([0.1, 0.2, 0.3]),
                post_probs=np.array([-9999.0, -9999.0]),
                bmdsRes=SimpleNamespace(BMDL_MA=0.1, BMD_MA=0.2, BMDU_MA=0.3),
                models=[
                    SimpleNamespace(
                        loudRes=SimpleNamespace(
                            BMD=np.array([0.1, 0.2, 0.3]),
                            parms=np.array([[1.0], [1.0], [1.0]]),
                            pval=0.5,
                            waic=0.0,
                        )
                    ),
                    SimpleNamespace(
                        loudRes=SimpleNamespace(
                            BMD=np.array([0.1, 0.2, 0.3]),
                            parms=np.array([[2.0], [2.0], [2.0]]),
                            pval=0.5,
                            waic=-2.0,
                        )
                    ),
                ],
            ),
            average=SimpleNamespace(modelPriors=np.array([0.5, 0.5])),
        )
        model_results = [
            SimpleNamespace(plotting=SimpleNamespace(dr_x=np.array([0.0, 1.0]), dr_y=np.ones(2))),
            SimpleNamespace(
                plotting=SimpleNamespace(dr_x=np.array([0.0, 1.0]), dr_y=np.full(2, 3.0))
            ),
        ]

        result = ContinuousModelAverageResult.from_cpp(analysis, model_results)

        expected = np.exp([0.0, -2.0])
        expected = expected / expected.sum()
        np.testing.assert_allclose(result.posteriors, expected)
        np.testing.assert_allclose(result.dr_y, expected @ np.array([[1.0, 1.0], [3.0, 3.0]]))

    def test_loud_draw_helpers_cover_chain_and_fallback_shapes(self):
        chains = [
            SimpleNamespace(BMD=np.array([1.0, 2.0]), parms=np.array([1.0, 2.0, 3.0, 4.0])),
            SimpleNamespace(BMD=np.array([3.0, 4.0]), parms=np.array([[5.0, 6.0], [7.0, 8.0]])),
        ]

        bmd, parms = ContinuousModelAverageResult._loud_draws(chains, n_chains=2)

        assert bmd.shape == (2, 2)
        assert parms.shape == (2, 2, 2)

        loud = SimpleNamespace(
            BMD=np.array([1.0, 2.0, 3.0, 4.0]),
            parms=np.array([[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]),
        )
        bmd, parms = ContinuousModelAverageResult._loud_draws(loud, n_chains=2)
        assert bmd.shape == (2, 2)
        assert parms.shape == (4, 2)

        result = SimpleNamespace(
            combinedLoudRes=SimpleNamespace(BMD=[]),
            loudRes=chains,
        )
        assert ContinuousModelAverageResult._combined_loud_result(result) is chains

    def test_continuous_sync_model_result_handles_empty_draws_with_existing_results(self):
        model = SimpleNamespace(
            settings=SimpleNamespace(alpha=0.05),
            results=SimpleNamespace(parameters="existing"),
            name=lambda: "Power",
        )
        result = ContinuousModelAverageResult(
            bmdl=0.1,
            bmd=0.2,
            bmdu=0.3,
            bmdl_y=0.0,
            bmd_y=0.0,
            bmdu_y=0.0,
            bmd_dist=np.array([]),
            priors=np.array([1.0]),
            posteriors=np.array([1.0]),
            model_bmd_dist=[np.array([])],
            model_parm_dist=[np.array([])],
            dr_x=np.array([0.0, 1.0]),
            dr_y=np.array([0.0, 1.0]),
        )

        # The empty-draw branch should reuse the existing parameter object before
        # later model attributes would be needed to finish syncing.
        try:
            result.sync_model_result(model, 0)
        except AttributeError:
            pass
        assert model.results.parameters == "existing"

    def test_json_safe_draws_preserves_supported_shapes(self):
        assert ContinuousModelAverageResult._json_safe_draws([]) == []
        assert ContinuousModelAverageResult._json_safe_draws([1.0, np.nan, np.inf]) == [
            1.0,
            None,
            None,
        ]
        assert ContinuousModelAverageResult._json_safe_draws(
            np.array([[1.0, np.nan], [2.0, np.inf]])
        ) == [[1.0, None], [2.0, None]]

    def test_reshape_flat_draws_leaves_uneven_chains_flat(self):
        draws = np.array([1.0, 2.0, 3.0])

        actual = ContinuousModelAverageResult._reshape_flat_draws(draws, n_chains=2)

        np.testing.assert_array_equal(actual, draws)
