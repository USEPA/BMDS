from dataclasses import dataclass
from unittest import mock

from pybmds.constants import BMDS_BLANK_VALUE, DistType, Dtype, LogicBin
from pybmds.recommender import checks
from pybmds.recommender.recommender import Rule, RuleClass


@dataclass
class ResultMock:
    pvalue: float
    dof: float

    def get_parameter(self, el: str):
        if el == "pvalue":
            return self.pvalue
        elif el == "dof":
            return self.dof


class TestChecks:
    def test_exists_rules(self, ddataset):
        model = mock.MagicMock()
        settings = Rule(rule_class=RuleClass.aic_missing, failure_bin=LogicBin.FAILURE)

        # good values
        for value in [-1, 0, 1]:
            model.results.get_parameter.return_value = value
            resp = checks.AicExists.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.NO_CHANGE
            assert resp.message == ""

        # bad values
        for value in [None, BMDS_BLANK_VALUE]:
            model.results.get_parameter.return_value = value
            resp = checks.AicExists.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.FAILURE
            assert resp.message == "AIC does not exist"

    def test_bmdl_exists_rules(self, ddataset):
        model = mock.MagicMock()
        settings = Rule(rule_class=RuleClass.bmdl_missing, failure_bin=LogicBin.FAILURE)

        # good values
        for value in [1e-8, 10]:
            model.results.get_parameter.return_value = value
            resp = checks.BmdlExists.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.NO_CHANGE
            assert resp.message == ""

        # special bad case for bmdl
        for value in [None, 0, BMDS_BLANK_VALUE]:
            model.results.get_parameter.return_value = value
            resp = checks.BmdlExists.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.FAILURE

    def test_greater_than_rules(self, ddataset):
        model = mock.MagicMock()
        settings = Rule(rule_class=RuleClass.gof, failure_bin=LogicBin.FAILURE, threshold=0.1)

        # good values
        for value in [BMDS_BLANK_VALUE, None, 0.1, 0.11]:
            model.results = ResultMock(value, 1)
            resp = checks.GoodnessOfFit.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.NO_CHANGE
            assert resp.message == ""

        # bad values
        for value in [0.09]:
            model.results = ResultMock(value, 1)
            resp = checks.GoodnessOfFit.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.FAILURE
            assert resp.message == "Goodness of fit p-value < 0.1"

    def test_less_than_rules(self, ddataset):
        dataset = mock.MagicMock()
        dataset.dtype = Dtype.DICHOTOMOUS
        model = mock.MagicMock()
        settings = Rule(rule_class=RuleClass.roi_large, failure_bin=LogicBin.FAILURE, threshold=2)

        # good values
        for value in [-2, 0, 2]:
            model.results.get_parameter.return_value = value
            resp = checks.LargeRoi.check(dataset, model, settings)
            assert resp.logic_bin == LogicBin.NO_CHANGE
            assert resp.message == ""

        # bad values
        for value in [-2.1, 2.1]:
            model.results.get_parameter.return_value = value
            resp = checks.LargeRoi.check(dataset, model, settings)
            assert resp.logic_bin == LogicBin.FAILURE
            assert resp.message == "|Residual near BMD| > 2.0"

    def test_gof(self, ddataset):
        model = mock.MagicMock()
        settings = Rule(rule_class=RuleClass.gof, failure_bin=LogicBin.FAILURE, threshold=0.1)

        # good values
        for p_value, df in [
            # these actually pass
            (0.1, 1),
            (0.11, 1),
            # skip test when we dont have a p_value
            (BMDS_BLANK_VALUE, 1),
            (None, 1),
            # skip test when DOF is 0ish
            (0.11, 0),
            (0.01, 0),
        ]:
            model.results = ResultMock(p_value, df)
            resp = checks.GoodnessOfFit.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.NO_CHANGE
            assert resp.message == ""

        # bad values
        for p_value, df in [(0.01, 1)]:
            model.results = ResultMock(p_value, df)
            resp = checks.GoodnessOfFit.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.FAILURE
            assert resp.message == "Goodness of fit p-value < 0.1"

    def test_variance_type(self, cdataset):
        model = mock.MagicMock()
        settings = Rule(
            rule_class=RuleClass.variance_type, failure_bin=LogicBin.FAILURE, threshold=0.1
        )
        for disttype, pvalue2, bin in [
            # if pvalue < 0.1 and disttype is normal, fail
            (DistType.normal, 0.09, LogicBin.FAILURE),
            (DistType.log_normal, 0.09, LogicBin.FAILURE),
            (DistType.normal, 0.11, LogicBin.NO_CHANGE),
            (DistType.log_normal, 0.11, LogicBin.NO_CHANGE),
            # if pvalue < 0.1 and disttype is nonconstant, pass
            (DistType.normal_ncv, 0.09, LogicBin.NO_CHANGE),
            (DistType.normal_ncv, 0.11, LogicBin.FAILURE),
        ]:
            model.settings.disttype = disttype
            model.results.tests.p_values = [None, pvalue2]
            resp = checks.VarianceType.check(cdataset, model, settings)
            assert resp.logic_bin == bin

    def test_variance_fit(self, cdataset):
        model = mock.MagicMock()
        settings = Rule(
            rule_class=RuleClass.variance_fit, failure_bin=LogicBin.FAILURE, threshold=0.1
        )
        # fmt: off
        for disttype, pvalues, bin, message in [
            # pass if test 2 > 0.1
            (DistType.normal, [1, 0.09, 1, 1], LogicBin.FAILURE, "Constant variance test failed (Test 2 p-value < 0.1)"),
            (DistType.log_normal, [1, 0.09, 1, 1], LogicBin.FAILURE, "Constant variance test failed (Test 2 p-value < 0.1)"),
            (DistType.normal, [0, 0.11, 0, 0], LogicBin.NO_CHANGE, ""),
            (DistType.log_normal, [0, 0.11, 0, 0], LogicBin.NO_CHANGE, ""),
            # pass if test 3 > 0.1
            (DistType.normal_ncv, [1, 1, 0.09, 1], LogicBin.FAILURE, "Nonconstant variance test failed (Test 3 p-value < 0.1)"),
            (DistType.normal_ncv, [0, 0, 0.11, 0], LogicBin.NO_CHANGE, ""),
        ]:
            # fmt: on
            model.settings.disttype = disttype
            model.results.tests.p_values = pvalues
            resp = checks.VarianceFit.check(cdataset, model, settings)
            assert resp.logic_bin == bin
            assert resp.message == message

    def test_zero_df(self, ddataset):
        model = mock.MagicMock()
        settings = Rule(rule_class=RuleClass.dof_zero, failure_bin=LogicBin.FAILURE)

        # good values
        for value in [0.1, 1]:
            model.results.get_parameter.return_value = value
            resp = checks.NoDegreesOfFreedom.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.NO_CHANGE
            assert resp.message == ""

        # bad values
        for value in [0]:
            model.results.get_parameter.return_value = value
            resp = checks.NoDegreesOfFreedom.check(ddataset, model, settings)
            assert resp.logic_bin == LogicBin.FAILURE
            assert resp.message == "Zero degrees of freedom; saturated model"
