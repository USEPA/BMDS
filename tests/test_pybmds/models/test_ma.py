import json

import numpy as np

import pybmds
from pybmds.constants import PriorClass


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
