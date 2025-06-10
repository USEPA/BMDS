// unit_tests.cpp
#include "unit_tests.h"

#include <vector>

#include "assert.h"
#include "bmds_helper.h"

int run_all_unitTests() {
  std::cout << "Running unit tests" << std::endl;
  objfunc_test();
  Nlogist_probs_test();
  // Nctr_probs_test();
  multitumor_ineq_constraint_test();
  multitumor_eq_constraint_test();
  dicho_AIC_penalty_test();
  cont_AIC_penalty_test();
  return 0;
}

void objfunc_test() {
  std::vector<double> x{1.5, 2.0, 3.2};
  std::vector<double> tmp;
  // assert(objfunc_bmdl(x, tmp, NULL)==1.5);
  expect_true(objfunc_bmdl(x, tmp, NULL) == 1.5);
}

void Nlogist_probs_test() {
  std::vector<double> Xi = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000};
  std::vector<double> Ls = {4,  5,  6,  6,  6,  6,  7,  7,  7,  9,  9,  9,  9,  9,  9,  11, 11, 11,
                            12, 12, 13, 13, 13, 14, 14, 2,  4,  5,  7,  7,  7,  7,  7,  8,  8,  9,
                            9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 14, 14, 3,  6,  6,  7,  7,
                            8,  8,  9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 13, 14, 14, 14, 14,
                            14, 14, 14, 6,  6,  6,  7,  7,  8,  8,  8,  10, 10, 10, 10, 10, 11, 11,
                            11, 12, 12, 13, 13, 13, 14, 14, 14, 3,  4,  4,  4,  5,  6,  6,  7,  7,
                            7,  7,  8,  9,  10, 10, 10, 11, 12, 12, 13, 13, 13, 14, 14, 14, 14};

  std::vector<double> p = {0.054599, -14.099836, 0.165150, -0.059984, 1.437543,
                           0.015216, 0.047678,   0.247516, 0.254625,  0.551903};
  // only care about Spec[0] and Spec[2]
  std::vector<bool> Spec = {false, false, false, false, false, false, false, false, false, false};
  bool compgrad = false;

  int Nobs = Xi.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));
  struct nestedObjData objData;
  objData.isBMDL = false;
  objData.smax = 14.0;
  objData.smin = 2.0;
  objData.Ls = Ls;
  objData.Xi = Xi;
  objData.Spec = Spec;
  objData.sijfixed = 9.28;
  objData.riskType = 1;
  objData.BMR = 0.1;
  objData.tD = 575.903266;

  std::vector<double> expProbs1 = {
      0.073024, 0.082237, 0.091449, 0.091449, 0.091449, 0.091449, 0.100662, 0.100662, 0.100662,
      0.119087, 0.119087, 0.119087, 0.119087, 0.119087, 0.119087, 0.137512, 0.137512, 0.137512,
      0.146725, 0.146725, 0.155937, 0.155937, 0.155937, 0.165150, 0.165150, 0.170946, 0.175634,
      0.178533, 0.185362, 0.185362, 0.185362, 0.185362, 0.185362, 0.189261, 0.189261, 0.193464,
      0.193464, 0.193464, 0.193464, 0.193464, 0.197956, 0.197956, 0.197956, 0.197956, 0.202723,
      0.207752, 0.213029, 0.218541, 0.218541, 0.310605, 0.300592, 0.300592, 0.298279, 0.298279,
      0.296467, 0.296467, 0.295146, 0.295146, 0.295146, 0.295146, 0.294306, 0.294306, 0.294306,
      0.294306, 0.293936, 0.294023, 0.294555, 0.294555, 0.295519, 0.295519, 0.295519, 0.295519,
      0.295519, 0.295519, 0.295519, 0.498023, 0.498023, 0.498023, 0.489821, 0.489821, 0.481997,
      0.481997, 0.481997, 0.467555, 0.467555, 0.467555, 0.467555, 0.467555, 0.460969, 0.460969,
      0.460969, 0.454823, 0.454823, 0.449128, 0.449128, 0.449128, 0.443894, 0.443894, 0.443894,
      0.797346, 0.789749, 0.789749, 0.789749, 0.782026, 0.774190, 0.774190, 0.766253, 0.766253,
      0.766253, 0.766253, 0.758230, 0.750136, 0.741988, 0.741988, 0.741988, 0.733804, 0.725603,
      0.725603, 0.717403, 0.717403, 0.717403, 0.709226, 0.709226, 0.709226, 0.709226
  };

  Nlogist_probs(probs, p, compgrad, gradij, &objData);
  for (int i = 0; i < Nobs; i++) {
    essentiallyEqual(expProbs1[i], probs[i], 1.5e-6);
  }
}

void Nctr_probs_test() {
  std::vector<double> Xi = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000};
  std::vector<double> Ls = {4,  5,  6,  6,  6,  6,  7,  7,  7,  9,  9,  9,  9,  9,  9,  11, 11, 11,
                            12, 12, 13, 13, 13, 14, 14, 2,  4,  5,  7,  7,  7,  7,  7,  8,  8,  9,
                            9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 14, 14, 3,  6,  6,  7,  7,
                            8,  8,  9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 13, 14, 14, 14, 14,
                            14, 14, 14, 6,  6,  6,  7,  7,  8,  8,  8,  10, 10, 10, 10, 10, 11, 11,
                            11, 12, 12, 13, 13, 13, 14, 14, 14, 3,  4,  4,  4,  5,  6,  6,  7,  7,
                            7,  7,  8,  9,  10, 10, 10, 11, 12, 12, 13, 13, 13, 14, 14, 14, 14};

  std::vector<double> p = {0.313819, 0.000030, 0.0,      -0.028465, 1.329996,
                           0.024342, 0.057298, 0.182997, 0.435914,  1.270733};

  // only care about Spec[0] and Spec[2]
  std::vector<bool> Spec = {false, false, false, false, false, false, false, false, false, false};
  bool compgrad = false;

  int Nobs = Xi.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));
  struct nestedObjData objData;
  objData.isBMDL = false;
  objData.smax = 4.464;
  objData.smin = -7.536;
  objData.smean = 9.536;
  objData.Ls = Ls;
  objData.Xi = Xi;
  objData.Spec = Spec;
  objData.sijfixed = 0;  // 9.28;
  objData.riskType = 1;
  objData.BMR = 0.1;
  objData.tD = 0.294277;

  std::vector<double> expProbs1 = {
      0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349,
      0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349,
      0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269349, 0.269351, 0.269351,
      0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350,
      0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350, 0.269350,
      0.269350, 0.269350, 0.269350, 0.269350, 0.269353, 0.269353, 0.269353, 0.269352, 0.269352,
      0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352,
      0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352, 0.269352,
      0.269352, 0.269352, 0.269352, 0.269358, 0.269358, 0.269358, 0.269358, 0.269358, 0.269358,
      0.269358, 0.269358, 0.269357, 0.269357, 0.269357, 0.269357, 0.269357, 0.269357, 0.269357,
      0.269357, 0.269357, 0.269357, 0.269356, 0.269356, 0.269356, 0.269356, 0.269356, 0.269356,
      0.269382, 0.269381, 0.269381, 0.269381, 0.269380, 0.269379, 0.269379, 0.269379, 0.269379,
      0.269379, 0.269379, 0.269378, 0.269377, 0.269376, 0.269376, 0.269376, 0.269376, 0.269375,
      0.269375, 0.269374, 0.269374, 0.269374, 0.269373, 0.269373, 0.269373, 0.269373
  };

  NCTR_probs(probs, p, compgrad, gradij, &objData);
  for (int i = 0; i < Nobs; i++) {
    essentiallyEqual(expProbs1[i], probs[i], 1.5e-6);
  }
}

void Nlogist_lk_test() {
  std::vector<double> Xi = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000};
  std::vector<double> Ls = {4,  5,  6,  6,  6,  6,  7,  7,  7,  9,  9,  9,  9,  9,  9,  11, 11, 11,
                            12, 12, 13, 13, 13, 14, 14, 2,  4,  5,  7,  7,  7,  7,  7,  8,  8,  9,
                            9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 14, 14, 3,  6,  6,  7,  7,
                            8,  8,  9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 13, 14, 14, 14, 14,
                            14, 14, 14, 6,  6,  6,  7,  7,  8,  8,  8,  10, 10, 10, 10, 10, 11, 11,
                            11, 12, 12, 13, 13, 13, 14, 14, 14, 3,  4,  4,  4,  5,  6,  6,  7,  7,
                            7,  7,  8,  9,  10, 10, 10, 11, 12, 12, 13, 13, 13, 14, 14, 14, 14};

  std::vector<double> p = {0.054599, -14.099836, 0.165150, -0.059984, 1.437543,
                           0.015216, 0.047678,   0.247516, 0.254625,  0.551903};
  // only care about Spec[0] and Spec[2]
  std::vector<bool> Spec = {false, false, false, false, false, false, false, false, false, false};
  bool compgrad = false;

  int Nobs = Xi.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));
  struct nestedObjData objData;
  objData.isBMDL = false;
  objData.smax = 14.0;
  objData.smin = 2.0;
  objData.Ls = Ls;
  objData.Xi = Xi;
  objData.Spec = Spec;
  objData.sijfixed = 9.28;
  objData.riskType = 1;
  objData.BMR = 0.1;
  objData.tD = 575.903266;

  double exp_lk = -337.540036;
  double lk = Nlogist_lk(p, &objData);

  essentiallyEqual(lk, exp_lk, 1.5e-6);
}

void multitumor_ineq_constraint_test() {
  const std::vector<double> x = {-2.8269246586783039, 1.3081031558978955E-021,
                                 0.69271407326031031, 1.1269537527262939E-021,
                                 -0.0000000000000000, 0.54665428796084770,
                                 0.16803414281323573, -0.0000000000000000,
                                 0.40754846974911302, 2.0785118843572268};

  double target = -103.61387669479440;
  std::vector<int> degree = {2, 2, 2};
  std::vector<double> doses1 = {0, 50, 100, 200, 400};
  std::vector<double> Y1 = {0, 1, 2, 10, 19};
  std::vector<double> n_group1 = {20, 20, 20, 20, 20};
  std::vector<double> doses2 = {0, 50, 100, 200, 400};
  std::vector<double> Y2 = {0, 1, 2, 4, 11};
  std::vector<double> n_group2 = {20, 20, 20, 20, 20};
  std::vector<double> doses3 = {0, 50, 100, 200, 400};
  std::vector<double> Y3 = {0, 2, 2, 6, 9};
  std::vector<double> n_group3 = {20, 20, 20, 20, 20};

  std::vector<std::vector<double>> doses;
  std::vector<std::vector<double>> Y;
  std::vector<std::vector<double>> n_group;
  doses.push_back(doses1);
  doses.push_back(doses2);
  doses.push_back(doses3);
  Y.push_back(Y1);
  Y.push_back(Y2);
  Y.push_back(Y3);
  n_group.push_back(n_group1);
  n_group.push_back(n_group2);
  n_group.push_back(n_group3);
  int nT = doses.size();

  // expected values come from BMDS-Model-Averaging repo results
  double expVal = 2.8421709430404E-013;
  const std::vector<double> expGrad = {
      0.0,
      -24.8065931556362,
      -3.08579588297736,
      -2.267848203526,
      -34.4556907297259,
      -3.08561840026068,
      -0.182575711582335,
      -37.6923046554479,
      -3.08575327396336,
      -0.18268740294875
  };

  double maxDose = 0;
  for (int i = 0; i < nT; i++) {
    double tmpMax = *std::max_element(doses[i].begin(), doses[i].end());
    if (tmpMax > maxDose) maxDose = tmpMax;
  }
  std::vector<double> grad(x.size());
  std::vector<int> nObs;
  struct msComboInEq ineq1;
  ineq1.nT = nT;
  ineq1.target = target;

  for (int i = 0; i < nT; i++) {
    std::vector<double> scaledDose = doses[i];
    nObs.push_back(scaledDose.size());
    for (int j = 0; j < scaledDose.size(); j++) {
      scaledDose[j] /= maxDose;
    }
    ineq1.doses.push_back(scaledDose);
    ineq1.Y.push_back(Y[i]);
    ineq1.n_group.push_back(n_group[i]);
  }
  ineq1.nObs = nObs;
  ineq1.degree = degree;

  double ineqVal = myInequalityConstraint1(x, grad, &ineq1);
  essentiallyEqual(ineqVal, expVal, 1e-18);

  for (int i = 0; i < grad.size(); i++) {
    essentiallyEqual(grad[i], expGrad[i], 1e-6);
  }
}

void multitumor_eq_constraint_test() {
  const std::vector<double> x = {-2.8269246586783039, 1.3081031558978955E-021,
                                 0.69271407326031031, 1.1269537527262939E-021,
                                 -0.0000000000000000, 0.54665428796084770,
                                 0.16803414281323573, -0.0000000000000000,
                                 0.40754846974911302, 2.0785118843572268};

  std::vector<int> degree = {2, 2, 2};
  int nT = degree.size();
  double bmr = 0.1;
  std::vector<double> grad(x.size());

  // expected values come from BMDS-Model-Averaging repo results
  double expVal = 6.38378239159465011e-16;
  const std::vector<double> expGrad = {0.113232419144430163,  0.0, 0.0591946176867892276,
                                       0.0035040027630851402, 0.0, 0.0591946176867892276,
                                       0.0035040027630851402, 0.0, 0.0591946176867892276,
                                       0.0035040027630851402};

  struct msComboEq eq1;
  eq1.bmr = bmr;
  eq1.nT = nT;
  eq1.degree = degree;

  double eqVal = myEqualityConstraint(x, grad, &eq1);

  essentiallyEqual(eqVal, expVal, 1e-18);

  for (int i = 0; i < grad.size(); i++) {
    essentiallyEqual(grad[i], expGrad[i], 1e-6);
  }
}

void cont_AIC_penalty_test (){
  struct continuous_analysis anal;
  anal.parms = 4;
  double prior[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1, 0.2, 1, -100, -100, 1, -18, 100, 100, 18,  18};
  anal.prior = prior;
  anal.model = power;

  struct BMDS_results bmdsRes;
  struct continuous_model_result res;
  double parm[4] = {5.05487, -0.0260177, 1, 0.698018};
  res.nparms = 4;
  res.parms = parm;
  res.max = 175.027;

  bool penalizeAIC;
  penalizeAIC = true;
  calcContAIC(&anal, &res, &bmdsRes, penalizeAIC);

  double AIC_penalized = bmdsRes.AIC;

  penalizeAIC = false;
  calcContAIC(&anal, &res, &bmdsRes, penalizeAIC);
  double AIC_unpenalized = bmdsRes.AIC;


  essentiallyEqual(AIC_unpenalized - 2,  AIC_penalized, 1e-6);
}

void dicho_AIC_penalty_test (){
  int estParmCount;
  struct dichotomous_analysis anal;
  anal.parms = 4;
  double prior[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 0, 0, 0, 18, 1e4, 1e4, 1e4};
  anal.prior = prior;
  anal.model = d_multistage;

  struct BMDS_results bmdsRes;
  struct dichotomous_model_result res;
  double parm[4] = {1.523e-8, 0, 1.05565e-5, 2.39079e-7};
  res.nparms = 4;
  res.parms = parm;
  res.max = 178.237;

  bool penalizeAIC;
  penalizeAIC = true;
  calcDichoAIC(&anal, &res, &bmdsRes, estParmCount, penalizeAIC);

  double AIC_penalized = bmdsRes.AIC;  

  penalizeAIC = false;
  calcDichoAIC(&anal, &res, &bmdsRes, estParmCount, penalizeAIC);
  double AIC_unpenalized = bmdsRes.AIC;


  essentiallyEqual(AIC_unpenalized - 2,  AIC_penalized, 1e-6);
}
