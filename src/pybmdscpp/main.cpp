
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bmds_helper.h"
#include "cmodeldefs.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_test1(py::module &);

////wrapper to expose vector of python_dichotomous_model_result structs to python
//struct python_dichotomousMA_resultPy : python_dichotomousMA_result{
//  py::list modelsPy;
//
//  public:
//  void modelsPy_set(const py::list & modelsPy)
//  //void my_setter(const py::list & models_vector_py)
//  {
//    std::cout<<"inside setter"<<std::endl;
//    this->models.clear();
//    for(auto& entry : modelsPy){
//      python_dichotomous_model_result res = entry.cast<python_dichotomous_model_result>;
//      this->models.push_back(res);
//    }
//    //for(int i = 0; i < py::len(modelsPy); ++i) {
//    //    //this->models.emplace_back(modelsPy[i].model);
//    //    //this->models.push_back(modelsPy[i].cast<python_dichotomous_model_result *>());
//    //}
//  }
//
//  py::list modelsPy_get()
//  {
//    py::list models_vector_py;
//      for(const auto & x : this->models) {
//          models_vector_py.append(x);
//      }
//    return models_vector_py;
//  }
//};


//pybind11::list python_dichotomousMA_resultPy::models_vector_py_get()
//{
//  pybind11::list models_vector_py;
//    for(const auto & x : this->models) {
//        models_vector_py.append(x);
//    }
//  return models_vector_py;
//}





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
       .def_readwrite("gof_chi_sqr_statistic", &python_dichotomous_model_result::gof_chi_sqr_statistic)
       .def_readwrite("gof", &python_dichotomous_model_result::gof)
       .def_readwrite("bmdsRes", &python_dichotomous_model_result::bmdsRes)
       .def_readwrite("aod", &python_dichotomous_model_result::aod);
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
    // functions
    init_test1(m);

}
