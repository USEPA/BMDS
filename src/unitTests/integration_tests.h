
#include "bmds_helper.h"

void runPythonDichoAnalysis();
void runPythonContAnalysis();
void runPythonMultitumorAnalysis();
int run_all_integrationTests();


std::vector<double> getMultitumorPrior(int degree, int prior_cols);
void printDichoModResult(struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes, bool showResultsOverride);
void createDichoAnalysisStructs(dich_model model, int modelType, bool restricted, int BMD_type, int degree, double BMR, double alpha, std::vector<double> &D, std::vector<double> &Y, std::vector<double> &N, python_dichotomous_analysis* anal, python_dichotomous_model_result* res);
