
#include "bmds_helper.h"

void runPythonDichoAnalysis();
void runPythonContAnalysis();
void runPythonMultitumorAnalysis();
int run_all_integrationTests();


std::vector<double> getMultitumorPrior(int degree, int prior_cols);
void printDichoModResult(struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes, bool showResultsOverride);
void printContModResult(struct python_continuous_analysis *pyAnal, struct python_continuous_model_result *pyRes, bool showResultsOverride);
void createDichoAnalysisStructs(dich_model model, int modelType, bool restricted, int BMD_type, int degree, double BMR, double alpha, std::vector<double> &D, std::vector<double> &Y, std::vector<double> &N, python_dichotomous_analysis* anal, python_dichotomous_model_result* res);
void createContAnalysisStructs(cont_model model, int modelType, bool restricted, distribution dist, bool detectAdvDir, bool isIncreasing, int BMD_type, double BMRF, int degree, double alpha, bool suffStat,  std::vector<double> &D, std::vector<double> &Y, std::vector<double> &N, std::vector<double> &SD, python_continuous_analysis* anal, python_continuous_model_result* res);
