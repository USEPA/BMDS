
#include "bmds_helper.h"

void runPythonDichoAnalysis();
void runPythonContAnalysis();
void runPythonMultitumorAnalysis();
int run_all_integrationTests();


std::vector<double> getMultitumorPrior(int degree, int prior_cols);
void printDichoModResult(struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes, bool showResultsOverride);
