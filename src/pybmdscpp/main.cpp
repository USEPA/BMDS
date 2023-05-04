
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bmds_helper.h"
#include "cmodeldefs.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_test1(py::module &);



PYBIND11_MODULE(bmdscore, m) {
    m.doc() = "bmdscore c++ interface";
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
    py::class_<test_struct>(m, "test_struct")
       .def(py::init<>())
       .def_readwrite("BMD", &test_struct::BMD)
       .def_readwrite("n", &test_struct::n)
       .def_readwrite("validResult", &test_struct::validResult)
       .def_readwrite("doses", &test_struct::doses);
    py::class_<python_dichotomous_analysis>(m, "python_dichotomous_analysis")
       .def(py::init<>())
       .def_readwrite("model", &python_dichotomous_analysis::model)
       .def_readwrite("n", &python_dichotomous_analysis::n)
       .def_readwrite("Y", &python_dichotomous_analysis::Y)
       .def_readwrite("doses", &python_dichotomous_analysis::doses)
       .def_readwrite("n_group", &python_dichotomous_analysis::n_group)
       .def_readwrite("prior", &python_dichotomous_analysis::prior)
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
       .def_readwrite("gof_chi_sqr_statistic", &python_dichotomous_model_result::gof_chi_sqr_statistic);
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
       .def_readwrite("validResult", &BMDS_results::validResult);
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
    // functions
    init_test1(m);

}
