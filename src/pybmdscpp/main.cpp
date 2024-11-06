
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bmds_helper.h"
#include "cmodeldefs.h"

namespace py = pybind11;

PYBIND11_MODULE(bmdscore, m) {
  m.doc() = "bmdscore C++ interface";

  m.def("version", &version, "bmdscore version number");

  py::enum_<dich_model>(m, "dich_model", py::arithmetic(), "Dichotomous model enumeration")
      .value("d_hill", d_hill, "dichotomous hill model enum")
      .value("d_gamma", d_gamma, "gamma model enum")
      .value("d_logistic", d_logistic, "logistic model enum")
      .value("d_loglogistic", d_loglogistic, "log-logistic model enum")
      .value("d_logprobit", d_logprobit, "log-probit model enum")
      .value("d_multistage", d_multistage, "multistage model enum")
      .value("d_probit", d_probit, "probit model enum")
      .value("d_qlinear", d_qlinear, "quantal linear model enum")
      .value("d_weibull", d_weibull, "weibull model enum")
      .export_values();

  py::enum_<nested_model>(m, "nested_model", py::arithmetic(), "Nested model enumeration")
      .value("nlogistic", nlogistic, "nested logistic model enum")
      .value("nctr", nctr, "NCTR model enum")
      .export_values();

  py::class_<test_struct>(m, "test_struct")
      .def(py::init<>())
      .def_readwrite("BMD", &test_struct::BMD)
      .def_readwrite("n", &test_struct::n)
      .def_readwrite("validResult", &test_struct::validResult)
      .def_readwrite("doses", &test_struct::doses);

  py::class_<dichotomous_GOF>(m, "dichotomous_GOF")
      .def(py::init<>())
      .def_readwrite("n", &dichotomous_GOF::n)
      .def_readwrite("expected", &dichotomous_GOF::expected)
      .def_readwrite("residual", &dichotomous_GOF::residual)
      .def_readwrite("test_statistic", &dichotomous_GOF::test_statistic)
      .def_readwrite("p_value", &dichotomous_GOF::p_value)
      .def_readwrite("df", &dichotomous_GOF::df)
      .def_readwrite("ebLower", &dichotomous_GOF::ebLower)
      .def_readwrite("ebUpper", &dichotomous_GOF::ebUpper);

  py::class_<BMDS_results>(m, "BMDS_results")
      .def(py::init<>())
      .def_readwrite("BMD", &BMDS_results::BMD)
      .def_readwrite("BMDL", &BMDS_results::BMDL)
      .def_readwrite("BMDU", &BMDS_results::BMDU)
      .def_readwrite("AIC", &BMDS_results::AIC)
      .def_readwrite("BIC_equiv", &BMDS_results::BIC_equiv)
      .def_readwrite("chisq", &BMDS_results::chisq)
      .def_readwrite("bounded", &BMDS_results::bounded)
      .def_readwrite("stdErr", &BMDS_results::stdErr)
      .def_readwrite("lowerConf", &BMDS_results::lowerConf)
      .def_readwrite("upperConf", &BMDS_results::upperConf)
      .def_readwrite("validResult", &BMDS_results::validResult)
      .def_readwrite("slopeFactor", &BMDS_results::slopeFactor)
      .def("setSlopeFactor", &BMDS_results::setSlopeFactor);

  py::class_<BMDSMA_results>(m, "BMDSMA_results")
      .def(py::init<>())
      .def_readwrite("BMD_MA", &BMDSMA_results::BMD_MA)
      .def_readwrite("BMDL_MA", &BMDSMA_results::BMDL_MA)
      .def_readwrite("BMDU_MA", &BMDSMA_results::BMDU_MA)
      .def_readwrite("BMD", &BMDSMA_results::BMD)
      .def_readwrite("BMDL", &BMDSMA_results::BMDL)
      .def_readwrite("BMDU", &BMDSMA_results::BMDU)
      .def_readwrite("ebLower", &BMDSMA_results::ebLower)
      .def_readwrite("ebUpper", &BMDSMA_results::ebUpper);

  py::class_<dicho_AOD>(m, "dicho_AOD")
      .def(py::init<>())
      .def_readwrite("fullLL", &dicho_AOD::fullLL)
      .def_readwrite("nFull", &dicho_AOD::nFull)
      .def_readwrite("redLL", &dicho_AOD::redLL)
      .def_readwrite("nRed", &dicho_AOD::nRed)
      .def_readwrite("fittedLL", &dicho_AOD::fittedLL)
      .def_readwrite("nFit", &dicho_AOD::nFit)
      .def_readwrite("devFit", &dicho_AOD::devFit)
      .def_readwrite("devRed", &dicho_AOD::devRed)
      .def_readwrite("dfFit", &dicho_AOD::dfFit)
      .def_readwrite("dfRed", &dicho_AOD::dfRed)
      .def_readwrite("pvFit", &dicho_AOD::pvFit)
      .def_readwrite("pvRed", &dicho_AOD::pvRed);

  py::enum_<cont_model>(m, "cont_model", py::arithmetic(), "Continuous model enumeration")
      .value("generic", generic, "generic model enum")
      .value("hill", hill, "continuous hill model enum")
      .value("exp_3", exp_3, "exponential 3 model enum")
      .value("exp_5", exp_5, "exponential 5 model enum")
      .value("power", power, "power model enum")
      .value("funl", funl, "funl model enum")
      .value("polynomial", polynomial, "polynomial model enum")
      .export_values();

  py::enum_<distribution>(m, "distribution", py::arithmetic(), "Continuous model distribution")
      .value("normal", normal, "normal distribution")
      .value("normal_ncv", normal_ncv, "normal non-constant variance distribution")
      .value("log_normal", log_normal, "log_normal distribution")
      .export_values();

  py::class_<continuous_GOF>(m, "continuous_GOF")
      .def(py::init<>())
      .def_readwrite("dose", &continuous_GOF::dose)
      .def_readwrite("size", &continuous_GOF::size)
      .def_readwrite("estMean", &continuous_GOF::estMean)
      .def_readwrite("calcMean", &continuous_GOF::calcMean)
      .def_readwrite("obsMean", &continuous_GOF::obsMean)
      .def_readwrite("estSD", &continuous_GOF::estSD)
      .def_readwrite("calcSD", &continuous_GOF::calcSD)
      .def_readwrite("obsSD", &continuous_GOF::obsSD)
      .def_readwrite("res", &continuous_GOF::res)
      .def_readwrite("n", &continuous_GOF::n)
      .def_readwrite("ebLower", &continuous_GOF::ebLower)
      .def_readwrite("ebUpper", &continuous_GOF::ebUpper);

  py::class_<continuous_AOD>(m, "continuous_AOD")
      .def(py::init<>())
      .def_readwrite("LL", &continuous_AOD::LL)
      .def_readwrite("nParms", &continuous_AOD::nParms)
      .def_readwrite("AIC", &continuous_AOD::AIC)
      .def_readwrite("addConst", &continuous_AOD::addConst)
      .def_readwrite("TOI", &continuous_AOD::TOI);

  py::class_<testsOfInterest>(m, "testsOfInterest")
      .def(py::init<>())
      .def_readwrite("llRatio", &testsOfInterest::llRatio)
      .def_readwrite("DF", &testsOfInterest::DF)
      .def_readwrite("pVal", &testsOfInterest::pVal);

  py::class_<nestedBootstrap>(m, "nestedBootstrap")
      .def(py::init<>())
      .def_readwrite("pVal", &nestedBootstrap::pVal)
      .def_readwrite("perc50", &nestedBootstrap::perc50)
      .def_readwrite("perc90", &nestedBootstrap::perc90)
      .def_readwrite("perc95", &nestedBootstrap::perc95)
      .def_readwrite("perc99", &nestedBootstrap::perc99);

  py::class_<nestedLitterData>(m, "nestedLitterData")
      .def(py::init<>())
      .def_readwrite("dose", &nestedLitterData::dose)
      .def_readwrite("LSC", &nestedLitterData::LSC)
      .def_readwrite("estProb", &nestedLitterData::estProb)
      .def_readwrite("litterSize", &nestedLitterData::litterSize)
      .def_readwrite("expected", &nestedLitterData::expected)
      .def_readwrite("observed", &nestedLitterData::observed)
      .def_readwrite("SR", &nestedLitterData::SR);

  py::class_<nestedReducedData>(m, "nestedReducedData")
      .def(py::init<>())
      .def_readwrite("dose", &nestedReducedData::dose)
      .def_readwrite("propAffect", &nestedReducedData::propAffect)
      .def_readwrite("lowerConf", &nestedReducedData::lowerConf)
      .def_readwrite("upperConf", &nestedReducedData::upperConf);

  py::class_<nestedSRData>(m, "nestedSRData")
      .def(py::init<>())
      .def_readwrite("minSR", &nestedSRData::minSR)
      .def_readwrite("avgSR", &nestedSRData::avgSR)
      .def_readwrite("maxSR", &nestedSRData::maxSR)
      .def_readwrite("minAbsSR", &nestedSRData::minAbsSR)
      .def_readwrite("avgAbsSR", &nestedSRData::avgAbsSR)
      .def_readwrite("maxAbsSR", &nestedSRData::maxAbsSR);

  py::class_<python_dichotomous_analysis>(m, "python_dichotomous_analysis")
      .def(py::init<>())
      .def_readwrite("model", &python_dichotomous_analysis::model)
      .def_readwrite("n", &python_dichotomous_analysis::n)
      .def_readwrite("Y", &python_dichotomous_analysis::Y)
      .def_readwrite("doses", &python_dichotomous_analysis::doses)
      .def_readwrite("n_group", &python_dichotomous_analysis::n_group)
      .def_readwrite("prior", &python_dichotomous_analysis::prior)
      .def_readwrite("BMR", &python_dichotomous_analysis::BMR)
      .def_readwrite("BMD_type", &python_dichotomous_analysis::BMD_type)
      .def_readwrite("alpha", &python_dichotomous_analysis::alpha)
      .def_readwrite("degree", &python_dichotomous_analysis::degree)
      .def_readwrite("samples", &python_dichotomous_analysis::samples)
      .def_readwrite("burnin", &python_dichotomous_analysis::burnin)
      .def_readwrite("parms", &python_dichotomous_analysis::parms)
      .def_readwrite("prior_cols", &python_dichotomous_analysis::prior_cols);

  py::class_<python_dichotomous_model_result>(m, "python_dichotomous_model_result")
      .def(py::init<>())
      .def_readwrite("model", &python_dichotomous_model_result::model)
      .def_readwrite("nparms", &python_dichotomous_model_result::nparms)
      .def_readwrite("parms", &python_dichotomous_model_result::parms)
      .def_readwrite("cov", &python_dichotomous_model_result::cov)
      .def_readwrite("max", &python_dichotomous_model_result::max)
      .def_readwrite("dist_numE", &python_dichotomous_model_result::dist_numE)
      .def_readwrite("model_df", &python_dichotomous_model_result::model_df)
      .def_readwrite("total_df", &python_dichotomous_model_result::total_df)
      .def_readwrite("bmd_dist", &python_dichotomous_model_result::bmd_dist)
      .def_readwrite("bmd", &python_dichotomous_model_result::bmd)
      .def_readwrite("gof_p_value", &python_dichotomous_model_result::gof_p_value)
      .def_readwrite("gof", &python_dichotomous_model_result::gof)
      .def_readwrite("bmdsRes", &python_dichotomous_model_result::bmdsRes)
      .def_readwrite("aod", &python_dichotomous_model_result::aod)
      .def_readwrite(
          "gof_chi_sqr_statistic", &python_dichotomous_model_result::gof_chi_sqr_statistic
      );

  py::class_<python_dichotomousMA_analysis>(m, "python_dichotomousMA_analysis")
      .def(py::init<>())
      .def_readwrite("nmodels", &python_dichotomousMA_analysis::nmodels)
      .def_readwrite("priors", &python_dichotomousMA_analysis::priors)
      .def_readwrite("nparms", &python_dichotomousMA_analysis::nparms)
      .def_readwrite("actual_parms", &python_dichotomousMA_analysis::actual_parms)
      .def_readwrite("prior_cols", &python_dichotomousMA_analysis::prior_cols)
      .def_readwrite("models", &python_dichotomousMA_analysis::models)
      .def_readwrite("modelPriors", &python_dichotomousMA_analysis::modelPriors)
      .def_readwrite("pyDA", &python_dichotomousMA_analysis::pyDA);

  py::class_<python_dichotomousMA_result>(m, "python_dichotomousMA_result")
      .def(py::init<>())
      .def_readwrite("nmodels", &python_dichotomousMA_result::nmodels)
      .def_readwrite("models", &python_dichotomousMA_result::models)
      .def_readwrite("dist_numE", &python_dichotomousMA_result::dist_numE)
      .def_readwrite("post_probs", &python_dichotomousMA_result::post_probs)
      .def_readwrite("bmd_dist", &python_dichotomousMA_result::bmd_dist)
      .def_readwrite("bmdsRes", &python_dichotomousMA_result::bmdsRes);

  py::class_<python_continuous_analysis>(m, "python_continuous_analysis")
      .def(py::init<>())
      .def_readwrite("model", &python_continuous_analysis::model)
      .def_readwrite("n", &python_continuous_analysis::n)
      .def_readwrite("suff_stat", &python_continuous_analysis::suff_stat)
      .def_readwrite("Y", &python_continuous_analysis::Y)
      .def_readwrite("doses", &python_continuous_analysis::doses)
      .def_readwrite("sd", &python_continuous_analysis::sd)
      .def_readwrite("n_group", &python_continuous_analysis::n_group)
      .def_readwrite("prior", &python_continuous_analysis::prior)
      .def_readwrite("BMD_type", &python_continuous_analysis::BMD_type)
      .def_readwrite("isIncreasing", &python_continuous_analysis::isIncreasing)
      .def_readwrite("BMR", &python_continuous_analysis::BMR)
      .def_readwrite("tail_prob", &python_continuous_analysis::tail_prob)
      .def_readwrite("disttype", &python_continuous_analysis::disttype)
      .def_readwrite("alpha", &python_continuous_analysis::alpha)
      .def_readwrite("samples", &python_continuous_analysis::samples)
      .def_readwrite("degree", &python_continuous_analysis::degree)
      .def_readwrite("burnin", &python_continuous_analysis::burnin)
      .def_readwrite("parms", &python_continuous_analysis::parms)
      .def_readwrite("prior_cols", &python_continuous_analysis::prior_cols)
      .def_readwrite("transform_dose", &python_continuous_analysis::transform_dose)
      .def_readwrite("restricted", &python_continuous_analysis::restricted)
      .def_readwrite("detectAdvDir", &python_continuous_analysis::detectAdvDir);

  py::class_<python_continuous_model_result>(m, "python_continuous_model_result")
      .def(py::init<>())
      .def_readwrite("model", &python_continuous_model_result::model)
      .def_readwrite("dist", &python_continuous_model_result::dist)
      .def_readwrite("nparms", &python_continuous_model_result::nparms)
      .def_readwrite("parms", &python_continuous_model_result::parms)
      .def_readwrite("cov", &python_continuous_model_result::cov)
      .def_readwrite("max", &python_continuous_model_result::max)
      .def_readwrite("dist_numE", &python_continuous_model_result::dist_numE)
      .def_readwrite("model_df", &python_continuous_model_result::model_df)
      .def_readwrite("total_df", &python_continuous_model_result::total_df)
      .def_readwrite("bmd", &python_continuous_model_result::bmd)
      .def_readwrite("bmd_dist", &python_continuous_model_result::bmd_dist)
      .def_readwrite("gof", &python_continuous_model_result::gof)
      .def_readwrite("bmdsRes", &python_continuous_model_result::bmdsRes)
      .def_readwrite("aod", &python_continuous_model_result::aod);

  py::class_<python_multitumor_analysis>(m, "python_multitumor_analysis")
      .def(py::init<>())
      .def_readwrite("ndatasets", &python_multitumor_analysis::ndatasets)
      .def_readwrite("models", &python_multitumor_analysis::models)
      .def_readwrite("n", &python_multitumor_analysis::n)
      .def_readwrite("nmodels", &python_multitumor_analysis::nmodels)
      .def_readwrite("BMD_type", &python_multitumor_analysis::BMD_type)
      .def_readwrite("BMR", &python_multitumor_analysis::BMR)
      .def_readwrite("alpha", &python_multitumor_analysis::alpha)
      .def_readwrite("prior_cols", &python_multitumor_analysis::prior_cols)
      .def_readwrite("degree", &python_multitumor_analysis::degree);

  py::class_<python_multitumor_result>(m, "python_multitumor_result")
      .def(py::init<>())
      .def_readwrite("ndatasets", &python_multitumor_result::ndatasets)
      .def_readwrite("validResult", &python_multitumor_result::validResult)
      .def_readwrite("nmodels", &python_multitumor_result::nmodels)
      .def_readwrite("models", &python_multitumor_result::models)
      .def_readwrite("selectedModelIndex", &python_multitumor_result::selectedModelIndex)
      .def_readwrite("BMD", &python_multitumor_result::BMD)
      .def_readwrite("BMDL", &python_multitumor_result::BMDL)
      .def_readwrite("BMDU", &python_multitumor_result::BMDU)
      .def_readwrite("slopeFactor", &python_multitumor_result::slopeFactor)
      .def_readwrite("combined_LL", &python_multitumor_result::combined_LL)
      .def_readwrite("combined_LL_const", &python_multitumor_result::combined_LL_const)
      .def("setSlopeFactor", &python_multitumor_result::setSlopeFactor);

  py::class_<python_nested_analysis>(m, "python_nested_analysis")
      .def(py::init<>())
      .def_readwrite("model", &python_nested_analysis::model)
      .def_readwrite("doses", &python_nested_analysis::doses)
      .def_readwrite("litterSize", &python_nested_analysis::litterSize)
      .def_readwrite("incidence", &python_nested_analysis::incidence)
      .def_readwrite("lsc", &python_nested_analysis::lsc)
      .def_readwrite("prior", &python_nested_analysis::prior)
      .def_readwrite("LSC_type", &python_nested_analysis::LSC_type)
      .def_readwrite("ILC_type", &python_nested_analysis::ILC_type)
      .def_readwrite("BMD_type", &python_nested_analysis::BMD_type)
      .def_readwrite("estBackground", &python_nested_analysis::estBackground)
      .def_readwrite("parms", &python_nested_analysis::parms)
      .def_readwrite("prior_cols", &python_nested_analysis::prior_cols)
      .def_readwrite("BMR", &python_nested_analysis::BMR)
      .def_readwrite("alpha", &python_nested_analysis::alpha)
      .def_readwrite("numBootRuns", &python_nested_analysis::numBootRuns)
      .def_readwrite("iterations", &python_nested_analysis::iterations)
      .def_readwrite("seed", &python_nested_analysis::seed);

  py::class_<python_nested_result>(m, "python_nested_result")
      .def(py::init<>())
      .def_readwrite("validResult", &python_nested_result::validResult)
      .def_readwrite("model", &python_nested_result::model)
      .def_readwrite("nparms", &python_nested_result::nparms)
      .def_readwrite("parms", &python_nested_result::parms)
      .def_readwrite("cov", &python_nested_result::cov)
      .def_readwrite("model_df", &python_nested_result::model_df)
      .def_readwrite("bmd", &python_nested_result::bmd)
      .def_readwrite("fixedLSC", &python_nested_result::fixedLSC)
      .def_readwrite("LL", &python_nested_result::LL)
      .def_readwrite("combPVal", &python_nested_result::combPVal)
      .def_readwrite("bmdsRes", &python_nested_result::bmdsRes)
      .def_readwrite("litter", &python_nested_result::litter)
      .def_readwrite("boot", &python_nested_result::boot)
      .def_readwrite("reduced", &python_nested_result::reduced)
      .def_readwrite("srData", &python_nested_result::srData);

  m.def(
      "pythonBMDSDicho", &pythonBMDSDicho, "Entry point to run BMDS dichotomous models",
      py::arg("python_dichotomous_analysis"), py::arg("python_dichotomous_model_result")
  );

  m.def(
      "pythonBMDSDichoMA", &pythonBMDSDichoMA, "Entry point to run BMDS dichotomous MA",
      py::arg("python_dichotomousMA_analysis"), py::arg("python_dichotomousMA_result")
  );

  m.def(
      "pythonBMDSCont", &pythonBMDSCont, "Entry point to run BMDS continuous models",
      py::arg("python_continuous_analysis"), py::arg("python_continuous_model_result")
  );

  m.def(
      "pythonBMDSMultitumor", &pythonBMDSMultitumor, "Entry point to run Multitumor analysis",
      py::arg("python_multitumor_analysis"), py::arg("python_multitumor_result")
  );

  m.def(
      "pythonBMDSNested", &pythonBMDSNested, "Entry point to run Nested analysis",
      py::arg("python_nested_analysis"), py::arg("python_nested_result")
  );
}
