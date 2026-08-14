// unit_tests.cpp
#include "unit_tests.h"

#include <vector>

#include "assert.h"
#include "bmds_helper.h"
#include "functional_generalized.h"

int run_all_unitTests() {
  std::cout << "Running unit tests" << std::endl;
  objfunc_test();
  // Nlogist_probs_test();
  //// Nctr_probs_test();
  // multitumor_ineq_constraint_test();
  multitumor_eq_constraint_test();
  dicho_AIC_penalty_test();
  cont_AIC_penalty_test();
  nested_AIC_penalty_test();
  // cont_loud_model_fit_test();
  dicho_loud_model_fit_test();
  pivotal_pvalue_test();
  additional_cont_calcs_test();
  additional_dicho_calcs_test();
  // rg_dg_test();
  // bridge_sample_test();

  return 0;
}

void objfunc_test() {
  std::vector<double> x{1.5, 2.0, 3.2};
  std::vector<double> tmp;
  // assert(objfunc_bmdl(x, tmp, NULL)==1.5);
  expect_true(objfunc_bmdl(x, tmp, NULL) == 1.5);
}

void dicho_loud_model_fit_test() {
  Eigen::MatrixXd Y{{2, 10}, {0, 10}, {9, 10}, {10, 10}};

  Eigen::MatrixXd D{{0}, {0}, {0.25}, {0.5}, {1.0}};

  struct fitInput fitIn;
  fitIn.Y = Y;
  fitIn.doses = D;
  fitIn.datatype = loud_datatype::l_dichotomous;

  fitIn.bmdtype = 1;  // 1 = extra ; added otherwise
  fitIn.bmr = 0.1;

  int iter = 5;  // # of rows in R
  int numParms = -1;

  Eigen::VectorXd bmd(iter);
  Eigen::MatrixXd parms;

  struct fitResult fitOut;
  fitOut.BMD = bmd;

  // QLINEAR
  fitIn.model = dich_model::d_qlinear;
  const Eigen::MatrixXd R_qlinear{
      {1.4420926, 3.446933, 1.7256582},
      {1.6642870, 3.081530, 1.0784842},
      {1.7704071, 2.648559, 1.5812958},
      {0.8703697, 4.126367, 0.4432010},
      {0.8644686, 3.296346, 0.3468643}
  };

  const Eigen::MatrixXd parms_qlinear{
      {0.2180138, 1.0977652},
      {0.2857488, 1.3499620},
      {0.2950550, 0.9839231},
      {0.1599963, 2.3331505},
      {0.1917769, 2.3516868}
  };

  const Eigen::VectorXd BMD_qlinear{{0.09597727, 0.07804702, 0.10708207, 0.04515805, 0.04480210}};

  numParms = parms_qlinear.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_qlinear(&fitIn, &fitOut, R_qlinear);

  expect_true(parms_qlinear.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_qlinear.isApprox(fitOut.BMD, 1.5e-6));

  // Logistic
  fitIn.model = dich_model::d_logistic;
  const Eigen::MatrixXd R_logistic{
      {0.6726732, 1.679524, 0.4848529},
      {0.6941349, 2.382565, 0.4983413},
      {0.7223239, 2.342939, 0.2171483},
      {0.2957612, 2.701273, 0.3539663},
      {0.4041959, 2.978572, 0.2713297}

  };

  const Eigen::MatrixXd parms_logistic{
      {-1.168628, 2.747888},
      {-1.423194, 3.243522},
      {-1.265323, 3.912631},
      {-2.335061, 4.471238},
      {-2.084480, 4.607595}
  };

  const Eigen::VectorXd BMD_logistic{{0.1398612, 0.1395135, 0.1044659, 0.1822487, 0.1509225}};

  numParms = parms_logistic.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_logistic(&fitIn, &fitOut, R_logistic);

  expect_true(parms_logistic.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_logistic.isApprox(fitOut.BMD, 1.5e-6));

  // Probit
  fitIn.model = dich_model::d_probit;
  const Eigen::MatrixXd R_probit{
      {0.6942345, 3.4937914, 0.5898628},
      {0.5902350, 3.5207400, 0.6307254},
      {0.8536257, 4.6657479, 1.0910242},
      {1.2049330, 5.2710791, 0.8344146},
      {0.9692136, 0.9741358, 0.6661034}
  };

  const Eigen::MatrixXd parms_probit{
      {-1.0567996, 2.214678},
      {-1.1528913, 2.265135},
      {-1.1304953, 2.104421},
      {-0.9748236, 2.179623},
      {-0.3280840, 0.986095}
  };

  const Eigen::VectorXd BMD_probit{{0.1447121, 0.1560583, 0.1641656, 0.1353929, 0.1648932}};

  numParms = parms_probit.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_probit(&fitIn, &fitOut, R_probit);

  expect_true(parms_probit.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_probit.isApprox(fitOut.BMD, 1.5e-6));

  // Mstage2
  fitIn.model = dich_model::d_multistage;
  const Eigen::MatrixXd R_mstage2{
      {0.6788134, 3.604108, 1.0704448, 0.4912015},
      {0.9514867, 3.854400, 0.9993588, 0.2770512},
      {0.6459280, 5.656038, 0.6612409, 0.1375121},
      {0.6913127, 6.061237, 0.6292018, 0.1160686},
      {0.9342977, 5.710257, 0.5344535, 0.4637425}
  };

  const Eigen::MatrixXd parms_mstage2{
      {0.12680123, 0.7240601, 0.7499991},
      {0.16390121, 0.4378502, 1.1425445},
      {0.09276301, 0.3103545, 1.9465710},
      {0.09365158, 0.2743843, 2.0895984},
      {0.13014301, 1.1399928, 1.3182525}
  };

  const Eigen::VectorXd BMD_mstage2{{0.12842868, 0.16745746, 0.16621102, 0.16829386, 0.08421996}};

  numParms = parms_mstage2.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_mstage2(&fitIn, &fitOut, R_mstage2);

  expect_true(parms_mstage2.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_mstage2.isApprox(fitOut.BMD, 1.5e-6));

  // Loglogistic
  fitIn.model = dich_model::d_loglogistic;
  const Eigen::MatrixXd R_loglogistic{
      {0.4300198, 2.574260, 0.7936150, 1.627801},
      {0.8013843, 3.043759, 0.4048359, 2.250533},
      {0.7590375, 3.995049, 0.2609827, 2.294265},
      {1.4462175, 4.395465, 0.4198194, 2.461597},
      {1.2609728, 4.619500, 0.2187056, 2.356580}
  };

  const Eigen::MatrixXd parms_loglogistic{
      {0.1132258, 1.176719, 1.627801},
      {0.1885619, 2.017367, 2.250533},
      {0.1513513, 2.728357, 2.294265},
      {0.2309697, 2.348504, 2.461597},
      {0.2067447, 3.050315, 2.356580}
  };

  const Eigen::VectorXd BMD_loglogistic{{0.1258455, 0.1537068, 0.1168446, 0.1577635, 0.1078773}};

  numParms = parms_loglogistic.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_loglogistic(&fitIn, &fitOut, R_loglogistic);

  expect_true(parms_loglogistic.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_loglogistic.isApprox(fitOut.BMD, 1.5e-6));

  // Logprobit
  fitIn.model = dich_model::d_logprobit;
  const Eigen::MatrixXd R_logprobit{
      {0.2432275, 2.293367, 0.31586361, 1.663363},
      {0.3478518, 2.673354, 0.34672000, 1.717312},
      {0.3634574, 2.775057, 0.35165604, 1.602705},
      {0.3526258, 2.375655, 0.34214157, 1.566058},
      {0.4888883, 2.767868, 0.04354112, 2.509254}
  };

  const Eigen::MatrixXd parms_logprobit{
      {0.08526946, 1.169723, 1.663363},
      {0.10328370, 1.201363, 1.717312},
      {0.10413742, 1.213506, 1.602705},
      {0.11484603, 1.146040, 1.566058},
      {0.14813459, 2.157399, 2.509254}
  };

  const Eigen::VectorXd BMD_logprobit{{0.2290797, 0.2355542, 0.2108140, 0.2122206, 0.2539781}};

  numParms = parms_logprobit.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_logprobit(&fitIn, &fitOut, R_logprobit);

  expect_true(parms_logprobit.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_logprobit.isApprox(fitOut.BMD, 1.5e-6));

  // DHill
  fitIn.model = dich_model::d_hill;
  const Eigen::MatrixXd R_dhill{
      {0.1369812, 0.9424855, 0.9423848, 2.008455},
      {0.1561788, 0.9676349, 0.8850810, 1.667696},
      {0.1539500, 0.9671256, 0.8906318, 1.649271},
      {0.1684538, 0.9379079, 0.9077965, 1.629195},
      {0.2569901, 0.9323692, 0.4914366, 1.777714}
  };

  const Eigen::MatrixXd parms_dhill{
      {0.1369812, 0.9424855, 0.9423848, 2.008455},
      {0.1561788, 0.9676349, 0.8850810, 1.667696},
      {0.1539500, 0.9671256, 0.8906318, 1.649271},
      {0.1684538, 0.9379079, 0.9077965, 1.629195},
      {0.2569901, 0.9323692, 0.4914366, 1.777714}
  };

  const Eigen::VectorXd BMD_dhill{{0.2164673, 0.1610112, 0.1572855, 0.1553642, 0.2302745}};

  numParms = parms_dhill.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_dhill(&fitIn, &fitOut, R_dhill);

  expect_true(parms_dhill.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_dhill.isApprox(fitOut.BMD, 1.5e-6));

  // Weibull
  fitIn.model = dich_model::d_weibull;
  const Eigen::MatrixXd R_weibull{
      {0.6783885, 4.921979, 0.7402839, 2.262244},
      {0.6519385, 5.115077, 0.2039971, 2.200412},
      {0.6313350, 5.102036, 0.2686074, 2.084397},
      {0.1813064, 2.513915, 0.1258766, 2.475641},
      {0.3451823, 2.463386, 0.3229572, 2.252823}
  };

  const Eigen::MatrixXd parms_weibull{
      {0.10699035, 2.262244, 2.034545},
      {0.10918392, 2.200412, 3.260949},
      {0.10518781, 2.084397, 2.995452},
      {0.06426803, 2.475641, 3.043153},
      {0.11022816, 2.252823, 2.154966}
  };

  const Eigen::VectorXd BMD_weibull{{0.2701662, 0.2101606, 0.2006964, 0.2570353, 0.2619201}};

  numParms = 3;
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_weibull(&fitIn, &fitOut, R_weibull);

  expect_true(parms_weibull.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_weibull.isApprox(fitOut.BMD, 1.5e-6));

  // Dgamma
  fitIn.model = dich_model::d_gamma;
  const Eigen::MatrixXd R_dgamma{
      {2.525968, 4.259757, 1.1667072, 1.700972},
      {1.319188, 5.250598, 0.6553346, 2.335863},
      {1.346608, 3.702852, 0.7797684, 2.431674},
      {1.399347, 3.423509, 1.0257794, 2.131211},
      {1.679361, 3.545789, 0.9825722, 2.205226}
  };

  const Eigen::MatrixXd parms_dgamma{
      {0.3176346, 1.700972, 2.504417},
      {0.1825836, 2.335863, 4.243680},
      {0.2310097, 2.431674, 3.756709},
      {0.2392604, 2.131211, 2.970455},
      {0.2705277, 2.205226, 3.149770}
  };

  const Eigen::VectorXd BMD_dgamma{{0.1531364, 0.1678831, 0.2039658, 0.2022576, 0.2033906}};

  numParms = parms_dgamma.cols();
  parms.resize(iter, numParms);

  fitOut.parms = parms;
  fit_dgamma(&fitIn, &fitOut, R_dgamma);

  expect_true(parms_dgamma.isApprox(fitOut.parms, 1.5e-6));
  expect_true(BMD_dgamma.isApprox(fitOut.BMD, 1.5e-6));
}

// these tests expect alpha to be returned not ln(alpha) as in BMDS
void cont_loud_model_fit_test() {
  Eigen::MatrixXd Y{
      {10.61764, 0.8937421, 20},
      {11.54771, 1.0580638, 20},
      {12.20492, 1.3528275, 20},
      {14.73715, 1.0778448, 19},
      {15.85227, 0.8350199, 19}
  };

  Eigen::MatrixXd D{{0}, {0.125}, {0.25}, {0.5}, {1.0}};

  bool isIncreasing = true;

  struct fitInput fitIn;
  fitIn.Y = Y;
  fitIn.doses = D;
  fitIn.datatype == loud_datatype::l_summary;

  fitIn.sign = 1.0;

  // create fitIn for each dist type
  // normal CV
  fitIn.bmdtype = 2;  // stddev
  fitIn.bmr = 1.0;
  fitIn.dist = distribution::normal;
  const struct fitInput fitInCV_sd = fitIn;

  fitIn.bmdtype = 3;  // reldev
  fitIn.bmr = 0.1;
  fitIn.dist = distribution::normal;
  const struct fitInput fitInCV_rel = fitIn;

  // normal NCV
  fitIn.dist = distribution::normal_ncv;
  fitIn.bmdtype = 2;  // stddev
  fitIn.bmr = 1.0;
  const struct fitInput fitInNCV_sd = fitIn;

  fitIn.bmdtype = 3;  // reldev
  fitIn.bmr = 0.1;
  const struct fitInput fitInNCV_rel = fitIn;

  // Lognormal CV
  fitIn.bmdtype = 2;  // stdev
  fitIn.bmr = 1.0;
  fitIn.dist = distribution::log_normal;
  const struct fitInput fitInLogCV_sd = fitIn;

  fitIn.bmdtype = 3;  // reldev
  fitIn.bmr = 0.1;
  fitIn.dist = distribution::log_normal;
  const struct fitInput fitInLogCV_rel = fitIn;

  int iter = 5;  // # of rows in R
  int numParms = -1;
  Eigen::VectorXd bmd(iter);
  Eigen::MatrixXd parms;

  struct fitResult fitOut;
  fitOut.BMD = bmd;

  ////////////////////
  // POWER CV
  //////////////////////

  numParms = 4;
  parms.resize(iter, numParms);

  fitOut.parms = parms;

  const Eigen::MatrixXd R_powerCV{
      {10.14480, 17.58850, 1.257526, 0.5927195},
      {10.14479, 17.58852, 1.257526, 0.5929868},
      {10.14479, 17.58851, 1.257526, 0.5928710},
      {10.14479, 17.58851, 1.257526, 0.5928589},
      {10.14478, 17.58854, 1.257525, 0.5933468}
  };

  const Eigen::MatrixXd powerCVParms{
      {10.14480, 7.443703, 1.257526, 1.687139},
      {10.14479, 7.443726, 1.257526, 1.686378},
      {10.14479, 7.443716, 1.257526, 1.686707},
      {10.14479, 7.443715, 1.257526, 1.686742},
      {10.14478, 7.443757, 1.257525, 1.685355}
  };

  const Eigen::VectorXd powerCV_stdev_BMD{{0.2494936, 0.2494482, 0.2494678, 0.2494699, 0.2493870}};
  const Eigen::VectorXd powerCV_reldev_BMD{{0.2049782, 0.2049775, 0.2049778, 0.2049778, 0.2049766}};

  fit_cpower(&fitInCV_sd, &fitOut, R_powerCV);
  expect_true(powerCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(powerCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cpower(&fitInCV_rel, &fitOut, R_powerCV);
  expect_true(powerCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  //////////////////////
  // POWER NCV
  ////////////////////
  numParms = 5;
  const Eigen::MatrixXd R_powerNCV{
      {10.39001, 15.60302, 1.515874, 0.3839357, 0.1988224},
      {10.39067, 15.60438, 1.515864, 0.3808076, 0.1854814},
      {10.39022, 15.60346, 1.515871, 0.3826683, 0.1953143},
      {10.39023, 15.60347, 1.515872, 0.3827986, 0.1924049},
      {10.37943, 15.58127, 1.516049, 0.4276003, 0.3808668}
  };

  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd powerNCVParms{
      {10.39001, 5.213014, 1.515874, 1.6183734, exp(-2.8310805)},
      {10.39067, 5.213711, 1.515864, 1.7689691, exp(-3.1755334)},
      {10.39022, 5.213237, 1.515871, 1.6539914, exp(-2.9111840)},
      {10.39023, 5.213244, 1.515872, 1.6917375, exp(-2.9998844)},
      {10.37943, 5.201839, 1.516049, 0.2849009, exp(0.1829478)}
  };

  const Eigen::VectorXd powerNCV_stdev_BMD{{0.4613986, 0.4626020, 0.4618883, 0.4618363, 0.4459658}};
  const Eigen::VectorXd powerNCV_reldev_BMD{{0.3450709, 0.3450523, 0.3450649, 0.3450652, 0.3453699}
  };

  fit_cpower(&fitInNCV_sd, &fitOut, R_powerNCV);
  expect_true(powerNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(powerNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cpower(&fitInNCV_rel, &fitOut, R_powerNCV);
  expect_true(powerNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  //////////////////
  // EXP3 CV
  //////////////////

  numParms = 4;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_exp3CV{
      {10.32854, 16.40208, 0.7347232, 1.558972},
      {10.32846, 16.40198, 0.7346882, 1.558979},
      {10.32847, 16.40198, 0.7346895, 1.558979},
      {10.32850, 16.40203, 0.7347064, 1.558975},
      {10.32853, 16.40207, 0.7347186, 1.558973}

  };

  const Eigen::MatrixXd exp3CVParms{
      {10.32854, 0.3501024, 0.7347232, 0.6414485},
      {10.32846, 0.3500856, 0.7346882, 0.6414454},
      {10.32847, 0.3500863, 0.7346895, 0.6414455},
      {10.32850, 0.3500944, 0.7347064, 0.6414470},
      {10.32853, 0.3501001, 0.7347186, 0.6414481}
  };

  const Eigen::VectorXd exp3CV_stdev_BMD{
      {0.08359856, 0.08358896, 0.08358931, 0.08359396, 0.08359729}
  };
  const Eigen::VectorXd exp3CV_reldev_BMD{{0.1165080, 0.1164957, 0.1164962, 0.1165021, 0.1165064}};

  fit_cexp3(&fitInCV_sd, &fitOut, R_exp3CV);
  expect_true(exp3CVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(exp3CV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cexp3(&fitInCV_rel, &fitOut, R_exp3CV);
  expect_true(exp3CV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  //////////////////
  // EXP3 NCV
  //////////////////

  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_exp3NCV{
      {10.74145, 16.73801, 1.487489, 0.3137527, 1.984308},
      {10.74138, 16.73809, 1.487592, 0.3125614, 1.984146},
      {10.74145, 16.73800, 1.487479, 0.3138605, 1.984321},
      {10.74150, 16.73795, 1.487417, 0.3145841, 1.984420},
      {10.74081, 16.73880, 1.488431, 0.3027921, 1.982807}
  };

  const Eigen::MatrixXd exp3NCVParms{
      {10.74145, 0.5789798, 1.487489, -4.158105, exp(11.03095)},
      {10.74138, 0.5790119, 1.487592, -4.166388, exp(11.05439)},
      {10.74145, 0.5789768, 1.487479, -4.157356, exp(11.02883)},
      {10.74150, 0.5789575, 1.487417, -4.152342, exp(11.01464)},
      {10.74081, 0.5792735, 1.488431, -4.235545, exp(11.25011)}
  };

  const Eigen::VectorXd exp3NCV_stdev_BMD{{0.4905265, 0.4911257, 0.4904725, 0.4901101, 0.4961462}};
  const Eigen::VectorXd exp3NCV_reldev_BMD{{0.3556629, 0.3556820, 0.3556611, 0.3556496, 0.3558381}};

  fit_cexp3(&fitInNCV_sd, &fitOut, R_exp3NCV);
  expect_true(exp3NCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(exp3NCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cexp3(&fitInNCV_rel, &fitOut, R_exp3NCV);
  expect_true(exp3NCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////
  // EXP3 LOGCV
  ///////////////////

  numParms = 4;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_exp3LogCV{
      {2.074359, 1.075104, 2.342920, 1.049676},
      {2.074831, 1.087564, 2.408487, 1.019421},
      {2.074691, 1.083938, 2.393829, 1.028519},
      {2.074695, 1.083803, 2.395533, 1.027121},
      {2.073364, 1.049481, 2.213176, 1.066825}
  };

  const Eigen::MatrixXd exp3LogCVParms{
      {7.959443, 0.9996822, 2.342920, log(0.9526748)},
      {7.963203, 0.9946935, 2.408487, log(0.9809492)},
      {7.962085, 0.9961268, 2.393829, log(0.9722722)},
      {7.962117, 0.9961875, 2.395533, log(0.9735949)},
      {7.951526, 1.0107213, 2.213176, log(0.9373611)}
  };

  const Eigen::VectorXd exp3LogCV_stdev_BMD{{1.045573, 1.059276, 1.055103, 1.055457, 1.031126}};
  const Eigen::VectorXd exp3LogCV_reldev_BMD{{0.3667891, 0.3788358, 0.3760367, 0.3762765, 0.3420610}
  };

  fit_cexp3(&fitInLogCV_sd, &fitOut, R_exp3LogCV);
  expect_true(exp3LogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(exp3LogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << exp3LogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_cexp3(&fitInLogCV_rel, &fitOut, R_exp3LogCV);
  expect_true(exp3LogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // EXP5 CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_exp5CV{
      {11.30294, 16.44929, 1.676246, 1.844872, 0.8144760},
      {11.30151, 16.44486, 1.676092, 1.843879, 0.8130661},
      {11.30061, 16.44208, 1.675995, 1.843257, 0.8121842},
      {11.30428, 16.45343, 1.676391, 1.845799, 0.8157870},
      {11.30429, 16.45347, 1.676392, 1.845808, 0.8157998}
  };

  const Eigen::MatrixXd exp5CVParms{
      {11.30294, 1.676246, 1.492102, 1.844872, 1.227783},
      {11.30151, 1.676092, 1.491948, 1.843879, 1.229912},
      {11.30061, 1.675995, 1.491851, 1.843257, 1.231248},
      {11.30428, 1.676391, 1.492246, 1.845799, 1.225810},
      {11.30429, 1.676392, 1.492247, 1.845808, 1.225791}
  };

  const Eigen::VectorXd exp5CV_stdev_BMD{{0.2639523, 0.2640700, 0.2641440, 0.2638431, 0.2638420}};
  const Eigen::VectorXd exp5CV_reldev_BMD{{0.2671595, 0.2671194, 0.2670944, 0.2671968, 0.2671971}};

  fit_cexp5(&fitInCV_sd, &fitOut, R_exp5CV);
  expect_true(exp5CVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(exp5CV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cexp5(&fitInCV_rel, &fitOut, R_exp5CV);
  expect_true(exp5CV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // EXP5 NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_exp5NCV{
      {10.81164, 16.82500, 1.973901, 2.442058, 0.7836753, 0.6786467},
      {10.83474, 16.83954, 1.969843, 2.433451, 0.7368818, 0.6585834},
      {10.82757, 16.83544, 1.971251, 2.436323, 0.7606610, 0.6687888},
      {10.82954, 16.83859, 1.971680, 2.436824, 0.7560435, 0.6667998},
      {10.78528, 16.80840, 1.979022, 2.452975, 0.8463022, 0.7060934}
  };

  const Eigen::MatrixXd exp5NCVParms{
      {10.81164, 1.973901, 1.559089, 2.442058, 0.3253742, exp(-0.5308331)},
      {10.83474, 1.969843, 1.557273, 2.433451, 0.2547470, exp(-0.3016726)},
      {10.82757, 1.971251, 1.557869, 2.436323, 0.2916227, exp(-0.4211056)},
      {10.82954, 1.971680, 1.557864, 2.436824, 0.2845722, exp(-0.3982736)},
      {10.78528, 1.979022, 1.561160, 2.452975, 0.4082274, exp(-0.8039603)}
  };

  const Eigen::VectorXd exp5NCV_stdev_BMD{{0.2657438, 0.2695897, 0.2675978, 0.2679259, 0.2610929}};
  const Eigen::VectorXd exp5NCV_reldev_BMD{{0.2605094, 0.2608178, 0.2607098, 0.2606898, 0.2601718}};

  fit_cexp5(&fitInNCV_sd, &fitOut, R_exp5NCV);
  expect_true(exp5NCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(exp5NCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cexp5(&fitInNCV_rel, &fitOut, R_exp5NCV);
  expect_true(exp5NCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // EXP5 LOGCV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_exp5LogCV{
      {1.748923, 2.941923, 2.564370, 1.0683093, 0.7626099},
      {1.750912, 2.938902, 2.565635, 1.0621063, 0.7611598},
      {1.756245, 2.930798, 2.568815, 1.0435225, 0.7568228},
      {1.761533, 2.922763, 2.572834, 1.0251954, 0.7525469},
      {1.779133, 2.896016, 2.585508, 0.9598605, 0.7373478}
  };

  const Eigen::MatrixXd exp5LogCVParms{
      {5.748410, 2.564370, 3.456401, 1.0683093, log(1.311286)},
      {5.759851, 2.565635, 3.441258, 1.0621063, log(1.313785)},
      {5.790654, 2.568815, 3.401958, 1.0435225, log(1.321313)},
      {5.821355, 2.572834, 3.363387, 1.0251954, log(1.328821)},
      {5.924718, 2.585508, 3.241374, 0.9598605, log(1.356212)}
  };

  const Eigen::VectorXd exp5LogCV_stdev_BMD{{0, 0, 0, 0, 0}};
  const Eigen::VectorXd exp5LogCV_reldev_BMD{
      {0.01986262, 0.01962696, 0.01888421, 0.01815240, 0.01551494}
  };

  fit_cexp5(&fitInLogCV_sd, &fitOut, R_exp5LogCV);
  expect_true(exp5LogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(exp5LogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << exp5LogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_cexp5(&fitInLogCV_rel, &fitOut, R_exp5LogCV);
  expect_true(exp5LogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // HILL CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_hillCV{
      {10.65489, 15.80145, 1.303407, 1.207015, 0.6025620},
      {10.63469, 15.77304, 1.277068, 1.191348, 0.5850586},
      {10.66320, 15.81313, 1.321693, 1.217895, 0.6151066},
      {10.66071, 15.80963, 1.316052, 1.214539, 0.6116256},
      {10.68051, 15.83748, 1.351207, 1.235461, 0.6354614}
  };

  const Eigen::MatrixXd hillCVParms{
      {10.65489, 12.23286, 1.303407, 1.207015, 1.659580},
      {10.63469, 12.01476, 1.277068, 1.191348, 1.709231},
      {10.66320, 12.38305, 1.321693, 1.217895, 1.625735},
      {10.66071, 12.33642, 1.316052, 1.214539, 1.634987},
      {10.68051, 12.63687, 1.351207, 1.235461, 1.573660}
  };

  const Eigen::VectorXd hillCV_stdev_BMD{{0.2214349, 0.2185837, 0.2234682, 0.2227723, 0.2267159}};
  const Eigen::VectorXd hillCV_reldev_BMD{{0.1860748, 0.1803578, 0.1900336, 0.1888123, 0.1964371}};

  fit_chill(&fitInCV_sd, &fitOut, R_hillCV);
  expect_true(hillCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(hillCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_chill(&fitInCV_rel, &fitOut, R_hillCV);
  expect_true(hillCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // HILL NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_hillNCV{
      {12.09806, 16.27880, 2.428470, 2.827084, 0.4899766, 0.1933174},
      {12.09780, 16.27790, 2.434746, 2.823307, 0.5148909, 0.1958377},
      {12.09793, 16.27834, 2.431686, 2.825010, 0.4956233, 0.1938859},
      {12.09854, 16.28046, 2.416787, 2.832798, 0.4277507, 0.1871914},
      {12.09880, 16.28137, 2.410370, 2.835964, 0.4037174, 0.1848708}
  };

  const Eigen::MatrixXd hillNCVParms{
      {12.09806, 55.54020, 2.428470, 2.827084, 3.133312, exp(-7.098091)},
      {12.09780, 55.73431, 2.434746, 2.823307, 3.257137, exp(-7.456318)},
      {12.09793, 55.63348, 2.431686, 2.825010, 3.162207, exp(-7.181550)},
      {12.09854, 55.11642, 2.416787, 2.832798, 2.783641, exp(-6.090638)},
      {12.09880, 54.88379, 2.410370, 2.835964, 2.630579, exp(-5.651274)}
  };

  const Eigen::VectorXd hillNCV_stdev_BMD{{0.6714700, 0.6651754, 0.6699198, 0.6885916, 0.6960470}};
  const Eigen::VectorXd hillNCV_reldev_BMD{{0.6322214, 0.6319117, 0.6320472, 0.6326456, 0.6328802}};

  fit_chill(&fitInNCV_sd, &fitOut, R_hillNCV);
  expect_true(hillNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(hillNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_chill(&fitInNCV_rel, &fitOut, R_hillNCV);
  expect_true(hillNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // HILL_EFSA CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_hill_efsaCV{
      {10.02324, 16.36718, 1.309544, 0.9652845, 0.6401729},
      {10.02234, 16.36447, 1.314428, 0.9701303, 0.6409054},
      {10.02251, 16.36497, 1.313528, 0.9692338, 0.6407703},
      {10.02349, 16.36793, 1.308172, 0.9639273, 0.6399671},
      {10.02415, 16.36994, 1.304571, 0.9603604, 0.6394266}
  };

  const Eigen::MatrixXd hill_efsaCVParms{
      {10.02324, 1.309544, 2.454041, 0.9652845, 1.562078},
      {10.02234, 1.314428, 2.457803, 0.9701303, 1.560293},
      {10.02251, 1.313528, 2.457105, 0.9692338, 1.560622},
      {10.02349, 1.308172, 2.452990, 0.9639273, 1.562580},
      {10.02415, 1.304571, 2.450249, 0.9603604, 1.563901}
  };

  const Eigen::VectorXd hill_efsaCV_stdev_BMD{
      {0.1128135, 0.1142339, 0.1139710, 0.1124165, 0.1113742}
  };
  const Eigen::VectorXd hill_efsaCV_reldev_BMD{
      {0.08806238, 0.08933479, 0.08909910, 0.08770697, 0.08677508}
  };

  fit_chill_efsa(&fitInCV_sd, &fitOut, R_hill_efsaCV);
  expect_true(hill_efsaCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(hill_efsaCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_chill_efsa(&fitInCV_rel, &fitOut, R_hill_efsaCV);
  expect_true(hill_efsaCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // HILL_EFSA NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_hill_efsaNCV{
      {12.06272, 15.93768, 2.248943, 2.176747, 0.7693763, 1.345567},
      {12.06878, 15.94227, 2.249322, 2.174841, 0.7695717, 1.347485},
      {12.04847, 15.92769, 2.248243, 2.180897, 0.7689705, 1.341564},
      {12.01877, 15.90611, 2.246599, 2.189308, 0.7682755, 1.333770},
      {12.03548, 15.91830, 2.247538, 2.183485, 0.7687131, 1.338787}
  };

  const Eigen::MatrixXd hill_efsaNCVParms{
      {12.06272, 2.248943, 3.196183, 2.176747, -2.006671, 5.259027},
      {12.06878, 2.249322, 3.192045, 2.174841, -2.012419, 5.274096},
      {12.04847, 2.248243, 3.206254, 2.180897, -1.993892, 5.225376},
      {12.01877, 2.246599, 3.226221, 2.189308, -1.968419, 5.158022},
      {12.03548, 2.247538, 3.213355, 2.183485, -1.984189, 5.199418}
  };

  const Eigen::VectorXd hill_efsaNCV_stdev_BMD{
      {0.5409143, 0.5406568, 0.5414021, 0.5423950, 0.5416361}
  };
  const Eigen::VectorXd hill_efsaNCV_reldev_BMD{
      {0.5557919, 0.5557098, 0.5558760, 0.5560579, 0.5557653}
  };

  fit_chill_efsa(&fitInNCV_sd, &fitOut, R_hill_efsaNCV);
  expect_true(hill_efsaNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(hill_efsaNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_chill_efsa(&fitInNCV_rel, &fitOut, R_hill_efsaNCV);
  expect_true(hill_efsaNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // HILL_EFSA LOGCV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_hill_efsaLogCV{
      {3.276238, 2.649848, 4.119359, 0.9196827, 2.431257},
      {3.276094, 2.650190, 4.118264, 0.9215378, 2.427465},
      {3.276088, 2.650203, 4.118268, 0.9215637, 2.429198},
      {3.276111, 2.650148, 4.118428, 0.9212743, 2.429384},
      {3.276547, 2.649112, 4.121061, 0.9165915, 2.440786}
  };

  const Eigen::MatrixXd hill_efsaLogCVParms{
      {26.47598, 4.119359, -1.176881, 0.9196827, log(0.4113099)},
      {26.47217, 4.118264, -1.179743, 0.9215378, log(0.4119523)},
      {26.47202, 4.118268, -1.179761, 0.9215637, log(0.4116585)},
      {26.47263, 4.118428, -1.179314, 0.9212743, log(0.4116270)},
      {26.48417, 4.121061, -1.172655, 0.9165915, log(0.4097042)}
  };

  const Eigen::VectorXd hill_efsaLogCV_stdev_BMD{{0, 0, 0, 0, 0}};
  const Eigen::VectorXd hill_efsaLogCV_reldev_BMD{{0, 0, 0, 0, 0}};

  fit_chill_efsa(&fitInLogCV_sd, &fitOut, R_hill_efsaLogCV);
  expect_true(hill_efsaLogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(hill_efsaLogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << hill_efsaLogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_chill_efsa(&fitInLogCV_rel, &fitOut, R_hill_efsaLogCV);
  expect_true(hill_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << hill_efsaLogCV_reldev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  ////////////////////////
  // INVEXP_EFSA CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_invexp_efsaCV{
      {11.98911, 16.05690, 1.247049, 0.8404389, 0.7715304},
      {11.98980, 16.04211, 1.245858, 0.8401301, 0.7715945},
      {11.99212, 16.00228, 1.242667, 0.8393123, 0.7717808},
      {11.99827, 15.90499, 1.234895, 0.8373311, 0.7722501},
      {12.00184, 15.85106, 1.230587, 0.8362347, 0.7725127}
  };

  const Eigen::MatrixXd invexp_efsaCVParms{
      {11.98911, 1.247049, 2.180753, 0.8404389, 1.296125},
      {11.98980, 1.245858, 2.174788, 0.8401301, 1.296018},
      {11.99212, 1.242667, 2.158643, 0.8393123, 1.295705},
      {11.99827, 1.234895, 2.119443, 0.8373311, 1.294917},
      {12.00184, 1.230587, 2.097900, 0.8362347, 1.294477}
  };

  const Eigen::VectorXd invexp_efsaCV_stdev_BMD{
      {0.4328969, 0.4332877, 0.4343952, 0.4373076, 0.4390427}
  };
  const Eigen::VectorXd invexp_efsaCV_reldev_BMD{
      {0.4437100, 0.4441579, 0.4454322, 0.4487815, 0.4507741}
  };

  fit_cinvexp_efsa(&fitInCV_sd, &fitOut, R_invexp_efsaCV);
  expect_true(invexp_efsaCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(invexp_efsaCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cinvexp_efsa(&fitInCV_rel, &fitOut, R_invexp_efsaCV);
  expect_true(invexp_efsaCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // INVEXP_EFSA NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_invexp_efsaNCV{
      {11.22846, 15.75746, 1.923303, 0.9551436, 1.018443, 0.9784878},
      {11.26382, 15.61273, 1.904467, 0.9613752, 1.021073, 0.9676984},
      {11.27399, 15.57034, 1.899105, 0.9630831, 1.021791, 0.9648177},
      {11.31230, 15.40773, 1.878549, 0.9694984, 1.024481, 0.9532550},
      {11.31312, 15.40421, 1.878104, 0.9696368, 1.024539, 0.9530331}
  };

  const Eigen::MatrixXd invexp_efsaNCVParms{
      {11.22846, 1.923303, 3.760334, 0.9551436, 0.1181076, -0.3039127},
      {11.26382, 1.904467, 3.592954, 0.9613752, 0.1644426, -0.4190677},
      {11.27399, 1.899105, 3.545612, 0.9630831, 0.1776977, -0.4520295},
      {11.31230, 1.878549, 3.369144, 0.9694984, 0.2332171, -0.5899454},
      {11.31312, 1.878104, 3.365411, 0.9696368, 0.2343817, -0.5928441}
  };

  const Eigen::VectorXd invexp_efsaNCV_stdev_BMD{
      {0.5435531, 0.5497982, 0.5516785, 0.5591659, 0.5593337}
  };
  const Eigen::VectorXd invexp_efsaNCV_reldev_BMD{
      {0.5650121, 0.5725550, 0.5748240, 0.5838551, 0.5840574}
  };

  fit_cinvexp_efsa(&fitInNCV_sd, &fitOut, R_invexp_efsaNCV);
  expect_true(invexp_efsaNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(invexp_efsaNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cinvexp_efsa(&fitInNCV_rel, &fitOut, R_invexp_efsaNCV);
  expect_true(invexp_efsaNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // INVEXP_EFSA LOGCV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_invexp_efsaLogCV{
      {2.170858, 3.311444, 0.10712713, 4.113787, 2.593265},
      {2.170763, 3.312395, 0.08819036, 4.104144, 2.598430},
      {2.170594, 3.314081, 0.06188680, 4.091029, 2.605702},
      {2.170754, 3.312481, 0.09084669, 4.105468, 2.598325},
      {2.170643, 3.313593, 0.07261393, 4.096308, 2.603397}
  };

  const Eigen::MatrixXd invexp_efsaLogCVParms{
      {8.765806, 0.10712713, 3.369291, 4.113787, log(0.3856143)},
      {8.764970, 0.08819036, 3.328424, 4.104144, log(0.3848478)},
      {8.763489, 0.06188680, 3.274165, 4.091029, log(0.3837737)},
      {8.764894, 0.09084669, 3.334945, 4.105468, log(0.3848634)},
      {8.763918, 0.07261393, 3.296879, 4.096308, log(0.3841136)}
  };

  const Eigen::VectorXd invexp_efsaLogCV_stdev_BMD{
      {0.5745927, 0.5493159, 0.5053713, 0.5530470, 0.5247766}
  };
  const Eigen::VectorXd invexp_efsaLogCV_reldev_BMD{
      {0.4390802, 0.4185078, 0.3834326, 0.4215712, 0.3988745}
  };

  fit_cinvexp_efsa(&fitInLogCV_sd, &fitOut, R_invexp_efsaLogCV);
  expect_true(invexp_efsaLogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(invexp_efsaLogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << invexp_efsaLogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_cinvexp_efsa(&fitInLogCV_rel, &fitOut, R_invexp_efsaLogCV);
  expect_true(invexp_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // LOG_EFSA CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_log_efsaCV{
      {10.58525, 15.94881, 1.345632, 0.9285246, 1.310680},
      {10.58694, 15.95196, 1.368187, 0.9402430, 1.314292},
      {10.58702, 15.95210, 1.369171, 0.9407546, 1.314464},
      {10.58668, 15.95147, 1.364706, 0.9384331, 1.313758},
      {10.58711, 15.95227, 1.370397, 0.9413908, 1.314682}
  };

  const Eigen::MatrixXd log_efsaCVParms{
      {10.58525, 1.345632, 1.821613, 0.9285246, 0.7629625},
      {10.58694, 1.368187, 1.813357, 0.9402430, 0.7608660},
      {10.58702, 1.369171, 1.813005, 0.9407546, 0.7607663},
      {10.58668, 1.364706, 1.814609, 0.9384331, 0.7611754},
      {10.58711, 1.370397, 1.812567, 0.9413908, 0.7606402}
  };

  const Eigen::VectorXd log_efsaCV_stdev_BMD{{0.1831855, 0.1847818, 0.1848508, 0.1845345, 0.1849362}
  };
  const Eigen::VectorXd log_efsaCV_reldev_BMD{
      {0.2068047, 0.2085582, 0.2086346, 0.2082877, 0.2087294}
  };

  fit_clog_efsa(&fitInCV_sd, &fitOut, R_log_efsaCV);
  expect_true(log_efsaCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(log_efsaCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_clog_efsa(&fitInCV_rel, &fitOut, R_log_efsaCV);
  expect_true(log_efsaCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // INVEXP_EFSA NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_log_efsaNCV{
      {9.911227, 15.93157, 2.162918, 0.8056158, 0.8733970, 1.623454},
      {9.911244, 15.93160, 2.162817, 0.8062361, 0.8735866, 1.624643},
      {9.911298, 15.93168, 2.162492, 0.8082204, 0.8741931, 1.628446},
      {9.911319, 15.93171, 2.162371, 0.8089633, 0.8744201, 1.629870},
      {9.911319, 15.93171, 2.162369, 0.8089744, 0.8744235, 1.629898},
  };

  const Eigen::MatrixXd log_efsaNCVParms{
      {9.911227, 2.162918, 1.778969, 0.8056158, -1.306102, 3.131131},
      {9.911244, 2.162817, 1.778983, 0.8062361, -1.307188, 3.133405},
      {9.911298, 2.162492, 1.779027, 0.8082204, -1.310652, 3.140665},
      {9.911319, 2.162371, 1.779043, 0.8089633, -1.311947, 3.143378},
      {9.911319, 2.162369, 1.779043, 0.8089744, -1.311976, 3.143439}
  };

  const Eigen::VectorXd log_efsaNCV_stdev_BMD{
      {0.09961412, 0.09978695, 0.10034004, 0.10054720, 0.10055029}
  };
  const Eigen::VectorXd log_efsaNCV_reldev_BMD{
      {0.09391613, 0.09409143, 0.09465258, 0.09486281, 0.09486594}
  };

  fit_clog_efsa(&fitInNCV_sd, &fitOut, R_log_efsaNCV);
  expect_true(log_efsaNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(log_efsaNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_clog_efsa(&fitInNCV_rel, &fitOut, R_log_efsaNCV);
  expect_true(log_efsaNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // LOG_EFSA LOGCV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_log_efsaLogCV{
      {2.376082, 3.220251, 1.329088, 0.1593405, 0.07730185},
      {2.375996, 3.220705, 1.329261, 0.1598616, 0.07730838},
      {2.375866, 3.221389, 1.329521, 0.1609399, 0.07732189},
      {2.375347, 3.224360, 1.330652, 0.1634249, 0.07735304},
      {2.375352, 3.222875, 1.330088, 0.1626546, 0.07734338}
  };

  const Eigen::MatrixXd log_efsaLogCVParms{
      {10.76265, 1.329088, 3.166795, 0.1593405, log(12.93630)},
      {10.76173, 1.329261, 3.168671, 0.1598616, log(12.93521)},
      {10.76033, 1.329521, 3.171500, 0.1609399, log(12.93295)},
      {10.75474, 1.330652, 3.183645, 0.1634249, log(12.92774)},
      {10.75480, 1.330088, 3.178540, 0.1626546, log(12.92935)}
  };

  const Eigen::VectorXd log_efsaLogCV_stdev_BMD{{0, 0, 0, 0, 0}};
  const Eigen::VectorXd log_efsaLogCV_reldev_BMD{
      {4.329355e-06, 4.491934e-06, 4.853505e-06, 5.721835e-06, 5.455290e-06}
  };

  fit_clog_efsa(&fitInLogCV_sd, &fitOut, R_log_efsaLogCV);
  expect_true(log_efsaLogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(log_efsaLogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << log_efsaLogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_clog_efsa(&fitInLogCV_rel, &fitOut, R_log_efsaLogCV);
  expect_true(log_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-5));

  ////////////////////////
  // GAMMA_EFSA CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_gamma_efsaCV{
      {11.04337, 16.33731, 1.477950, 2.055071, 1.0166559},
      {10.99190, 16.25361, 1.475608, 2.041575, 0.9826547},
      {11.00410, 16.27322, 1.476393, 2.044560, 0.9859853},
      {11.02579, 16.30842, 1.477441, 2.050913, 1.0001329},
      {10.91135, 16.12829, 1.467345, 2.013359, 0.9077365}
  };

  const Eigen::MatrixXd gamma_efsaCVParms{
      {11.04337, 1.477950, 2.147723, 2.055071, 0.9836170},
      {10.99190, 1.475608, 2.136907, 2.041575, 1.0176515},
      {11.00410, 1.476393, 2.139021, 2.044560, 1.0142140},
      {11.02579, 1.477441, 2.144063, 2.050913, 0.9998671},
      {10.91135, 1.467345, 2.119688, 2.013359, 1.1016413}
  };

  const Eigen::VectorXd gamma_efsaCV_stdev_BMD{
      {0.3288132, 0.3308554, 0.3307659, 0.3300321, 0.3354937}
  };
  const Eigen::VectorXd gamma_efsaCV_reldev_BMD{
      {0.3498703, 0.3477864, 0.3482332, 0.3492353, 0.3432667}
  };

  fit_cgamma_efsa(&fitInCV_sd, &fitOut, R_gamma_efsaCV);
  expect_true(gamma_efsaCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(gamma_efsaCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cgamma_efsa(&fitInCV_rel, &fitOut, R_gamma_efsaCV);
  expect_true(gamma_efsaCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // GAMMA_EFSA NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_gamma_efsaNCV{
      {10.58641, 16.01285, 0.8623874, 1.0495368, 1.470493, 1.030325},
      {10.58404, 16.01223, 0.8703607, 1.0256239, 1.480738, 1.025680},
      {10.58325, 16.01203, 0.8729801, 1.0186506, 1.483431, 1.024364},
      {10.58137, 16.01154, 0.8793035, 1.0001101, 1.491611, 1.020737},
      {10.57975, 16.01112, 0.8848159, 0.9836602, 1.498777, 1.017495}

  };

  const Eigen::MatrixXd gamma_efsaNCVParms{
      {10.58641, 0.8623874, 1.923162, 1.0495368, 0.8596095, -2.413907},
      {10.58404, 0.8703607, 1.900538, 1.0256239, 0.8869089, -2.485067},
      {10.58325, 0.8729801, 1.893987, 1.0186506, 0.8942651, -2.504173},
      {10.58137, 0.8793035, 1.877417, 1.0001101, 0.9157869, -2.560285},
      {10.57975, 0.8848159, 1.863314, 0.9836602, 0.9347479, -2.609665}
  };

  const Eigen::VectorXd gamma_efsaNCV_stdev_BMD{
      {0.1179187, 0.1116466, 0.1098306, 0.1050318, 0.1008303}
  };
  const Eigen::VectorXd gamma_efsaNCV_reldev_BMD{
      {0.1517057, 0.1449052, 0.1429083, 0.1376565, 0.1330235}
  };

  fit_cgamma_efsa(&fitInNCV_sd, &fitOut, R_gamma_efsaNCV);
  expect_true(gamma_efsaNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(gamma_efsaNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_cgamma_efsa(&fitInNCV_rel, &fitOut, R_gamma_efsaNCV);
  expect_true(gamma_efsaNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // GAMMA_EFSA LOGCV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_gamma_efsaLogCV{
      {3.412436, 1.931381, 1.953997, 0.5720145, 2.949811},
      {3.412436, 1.931362, 1.953983, 0.5719660, 2.949804},
      {3.412437, 1.931235, 1.953895, 0.5716567, 2.949756},
      {3.412437, 1.931264, 1.953915, 0.5717267, 2.949767},
      {3.412437, 1.931245, 1.953901, 0.5716806, 2.949760}
  };

  const Eigen::MatrixXd gamma_efsaLogCVParms{
      {30.33905, 1.953997, 0.1791335, 0.5720145, log(0.3390048)},
      {30.33906, 1.953983, 0.1791344, 0.5719660, log(0.3390056)},
      {30.33910, 1.953895, 0.1791402, 0.5716567, log(0.3390111)},
      {30.33909, 1.953915, 0.1791389, 0.5717267, log(0.3390099)},
      {30.33910, 1.953901, 0.1791398, 0.5716806, log(0.3390107)}
  };

  const Eigen::VectorXd gamma_efsaLogCV_stdev_BMD{{0, 0, 0, 0, 0}};
  const Eigen::VectorXd gamma_efsaLogCV_reldev_BMD{{0, 0, 0, 0, 0}};

  fit_cgamma_efsa(&fitInLogCV_sd, &fitOut, R_gamma_efsaLogCV);
  expect_true(gamma_efsaLogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(gamma_efsaLogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << gamma_efsaLogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_cgamma_efsa(&fitInLogCV_rel, &fitOut, R_gamma_efsaLogCV);
  expect_true(gamma_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << gamma_efsaLogCV_reldev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  ////////////////////////
  // LMS_EFSA CV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_lms_efsaCV{
      {10.76268, 15.64665, 1.154099, 2.067170, 0.8248621},
      {10.76287, 15.64746, 1.154211, 2.067450, 0.8180637},
      {10.76015, 15.63588, 1.152598, 2.063401, 0.9102798},
      {10.76025, 15.63633, 1.152661, 2.063560, 0.9061703},
      {10.76037, 15.63682, 1.152729, 2.063732, 0.9045700}
  };

  const Eigen::MatrixXd lms_efsaCVParms{
      {10.76268, 1.154099, 1.472649, 2.067170, 1.212324},
      {10.76287, 1.154211, 1.472692, 2.067450, 1.222399},
      {10.76015, 1.152598, 1.472066, 2.063401, 1.098563},
      {10.76025, 1.152661, 1.472090, 2.063560, 1.103545},
      {10.76037, 1.152729, 1.472117, 2.063732, 1.105498}
  };

  const Eigen::VectorXd lms_efsaCV_stdev_BMD{{0.1634784, 0.1640736, 0.1566035, 0.1569142, 0.1570243}
  };
  const Eigen::VectorXd lms_efsaCV_reldev_BMD{
      {0.1600775, 0.1600502, 0.1604460, 0.1604305, 0.1604137}
  };

  fit_clms_efsa(&fitInCV_sd, &fitOut, R_lms_efsaCV);
  expect_true(lms_efsaCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(lms_efsaCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_clms_efsa(&fitInCV_rel, &fitOut, R_lms_efsaCV);
  expect_true(lms_efsaCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // LMS_EFSA NCV
  ///////////////////////
  numParms = 6;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_lms_efsaNCV{
      {9.732914, 15.64322, 2.287442, 2.173154, 1.214485, 1.102238},
      {9.734242, 15.64531, 2.293258, 2.181287, 1.200332, 1.088552},
      {9.733977, 15.64489, 2.292082, 2.179651, 1.203153, 1.091283},
      {9.733981, 15.64490, 2.292100, 2.179688, 1.203160, 1.091279},
      {9.733304, 15.64383, 2.289146, 2.175557, 1.209988, 1.097878}
  };

  const Eigen::MatrixXd lms_efsaNCVParms{
      {9.732914, 2.287442, 1.614348, 2.173154, 0.2043692, -0.6593654},
      {9.734242, 2.293258, 1.614244, 2.181287, 0.2059958, -0.6513722},
      {9.733977, 2.292082, 1.614265, 2.179651, 0.2056643, -0.6529601},
      {9.733981, 2.292100, 1.614265, 2.179688, 0.2056826, -0.6530075},
      {9.733304, 2.289146, 1.614317, 2.175557, 0.2049024, -0.6568769}
  };

  const Eigen::VectorXd lms_efsaNCV_stdev_BMD{
      {0.06760979, 0.06785231, 0.06780391, 0.06780321, 0.06769025}
  };
  const Eigen::VectorXd lms_efsaNCV_reldev_BMD{
      {0.07265337, 0.07248765, 0.07252109, 0.07252057, 0.07260473}
  };

  fit_clms_efsa(&fitInNCV_sd, &fitOut, R_lms_efsaNCV);
  expect_true(lms_efsaNCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(lms_efsaNCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  fit_clms_efsa(&fitInNCV_rel, &fitOut, R_lms_efsaNCV);
  expect_true(lms_efsaNCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));

  ////////////////////////
  // LMS_EFSA LOGCV
  ///////////////////////
  numParms = 5;
  parms.resize(iter, numParms);
  fitOut.parms = parms;

  const Eigen::MatrixXd R_lms_efsaLogCV{
      {2.565590, 2.720145, 1.294782, 2.019590, 0.3270511},
      {2.565251, 2.719362, 1.299844, 2.021746, 0.3279970},
      {2.569078, 2.728208, 1.244163, 1.999033, 0.3174488},
      {2.568165, 2.726098, 1.257789, 2.004633, 0.3200272},
      {2.564318, 2.717209, 1.313767, 2.027419, 0.3306319}
  };

  const Eigen::MatrixXd lms_efsaLogCVParms{
      {13.00834, 1.294782, 1.173444, 2.019590, log(3.057626)},
      {13.00392, 1.299844, 1.172860, 2.021746, log(3.048809)},
      {13.05378, 1.244163, 1.179499, 1.999033, log(3.150114)},
      {13.04187, 1.257789, 1.177901, 2.004633, log(3.124734)},
      {12.99180, 1.313767, 1.171260, 2.027419, log(3.024512)}
  };

  const Eigen::VectorXd lms_efsaLogCV_stdev_BMD{{0, 0, 0, 0, 0}};
  const Eigen::VectorXd lms_efsaLogCV_reldev_BMD{
      {0.4062564, 0.4070011, 0.3989147, 0.4006745, 0.4091239}
  };

  fit_clms_efsa(&fitInLogCV_sd, &fitOut, R_lms_efsaLogCV);
  expect_true(lms_efsaLogCVParms.isApprox(fitOut.parms, 1.5e-6));
  expect_true(lms_efsaLogCV_stdev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << lms_efsaLogCV_stdev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  fit_clms_efsa(&fitInLogCV_rel, &fitOut, R_lms_efsaLogCV);
  expect_true(lms_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-5));
}

void pivotal_pvalue_test() {
  // Individual lognormal diagnostics use a model prediction on the log scale.
  // Verify that both the pivotal statistic and pointwise likelihood use that
  // scale and standardize residuals by the log-scale variance, not its SD.
  Eigen::MatrixXd Y_log{{std::exp(2.0)}, {std::exp(2.2)}};
  Eigen::VectorXd mu_log{{2.1}, {2.1}};
  Eigen::VectorXd parms_log{{1.0}, {100.0}};  // final parameter is precision
  double q_log =
      getQVals(Y_log, parms_log, mu_log, distribution::log_normal, loud_datatype::l_individual);
  expect_true(essentiallyEqual(2.0, q_log, 1e-12));

  Eigen::VectorXd loglik_log = loud_likelihood(
      Y_log, parms_log, mu_log, distribution::log_normal, loud_datatype::l_individual
  );
  double pi = std::acos(-1.0);
  double expected_log_density_0 = -2.0 - 0.5 * std::log(2.0 * pi * 0.01) - 0.5;
  double expected_log_density_1 = -2.2 - 0.5 * std::log(2.0 * pi * 0.01) - 0.5;
  expect_true(essentiallyEqual(expected_log_density_0, loglik_log(0), 1e-12));
  expect_true(essentiallyEqual(expected_log_density_1, loglik_log(1), 1e-12));

  // individual data
  Eigen::MatrixXd Y{{11.7692}, {11.5889}, {10.9485}, {12.0335}, {10.5843}, {10.9571}, {12.1292},
                    {12.5896}, {12.9235}, {11.1651}, {9.2031},  {11.7029}, {11.9342}, {11.8759},
                    {10.6886}, {12.0227}, {11.277},  {10.811},  {14.2668}, {10.4833}, {11.4515},
                    {12.5341}, {12.4558}, {13.7429}, {13.0073}, {12.5508}, {10.4007}, {12.8715},
                    {13.4137}, {9.33564}, {10.9298}, {10.6283}, {12.9618}, {12.2369}, {14.8123},
                    {11.0151}, {13.5539}, {10.7161}, {12.7046}, {12.7758}, {14.4463}, {14.8051},
                    {17.1731}, {14.4799}, {14.8711}, {13.9808}, {12.8826}, {15.3315}, {14.9358},
                    {14.0363}, {14.8257}, {14.0969}, {15.5807}, {15.202},  {14.1877}, {14.5575},
                    {13.1242}, {17.0857}, {14.4028}};

  //  Eigen::MatrixXd D{{0}, {0.125}, {0.25}, {0.5}, {1.0}};
  Eigen::MatrixXd D{{0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125},
                    {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125},
                    {0.125}, {0.125}, {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},
                    {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},
                    {0.25},  {0.25},  {0.25},  {0.25},  {0.5},   {0.5},   {0.5},   {0.5},   {0.5},
                    {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},
                    {0.5},   {0.5},   {0.5},   {0.5},   {0.5}};

  // Individual CV
  struct fitInput loudIn;
  loudIn.model = cont_model::power;
  loudIn.dist = distribution::normal;
  loudIn.datatype = loud_datatype::l_individual;
  loudIn.iter = 5;
  loudIn.burnin = 2;
  loudIn.qlev = 0.9;
  loudIn.df_override = BMDS_MISSING;
  loudIn.Y = Y;
  loudIn.doses = D;

  Eigen::MatrixXd R_power_cv{
      {10.14480, 17.58850, 1.257526, 0.5927195},
      {10.14479, 17.58852, 1.257526, 0.5929868},
      {10.14479, 17.58851, 1.257526, 0.5928710},
      {10.14479, 17.58851, 1.257526, 0.5928589},
      {10.14478, 17.58854, 1.257525, 0.5933468}
  };

  int S = R_power_cv.rows();
  int N = loudIn.Y.rows();
  Eigen::MatrixXd mu_power_cv(S, N);
  int model_typ = getLoudModelType(loudIn.model, loudIn.dist, loudIn.datatype);
  ptr2 model_fun = choose_nonlinearity2(model_typ);
  for (int i = 0; i < S; i++) {
    mu_power_cv.row(i) = model_fun(R_power_cv.row(i), loudIn.doses);
  }

  double expected_pval = 0.05070569;
  double returned_pval = pivotal_pvalue(R_power_cv, &loudIn, mu_power_cv);
  expect_true(essentiallyEqual(expected_pval, returned_pval, 0.001));

  // R contains retained post-burn-in draws, so the configured warm-up length
  // must not trim the draw matrix again or change the pivotal p-value.
  loudIn.burnin = 5;
  double returned_pval_different_burnin = pivotal_pvalue(R_power_cv, &loudIn, mu_power_cv);
  expect_true(essentiallyEqual(returned_pval, returned_pval_different_burnin, 1e-12));
  loudIn.burnin = 2;

  // Individual NCV
  loudIn.dist = distribution::normal_ncv;
  Eigen::MatrixXd R_power_ncv{
      {10.39001, 15.60302, 1.515874, 0.3839357, 0.1988224},
      {10.39067, 15.60438, 1.515864, 0.3808076, 0.1854814},
      {10.39022, 15.60346, 1.515871, 0.3826683, 0.1953143},
      {10.39023, 15.60347, 1.515872, 0.3827986, 0.1924049},
      {10.37943, 15.58127, 1.516049, 0.4276003, 0.3808668}
  };

  // int S = R_power_cv.rows();
  // int N = loudIn.Y.rows();
  Eigen::MatrixXd mu_power_ncv(S, N);
  model_typ = getLoudModelType(loudIn.model, loudIn.dist, loudIn.datatype);
  model_fun = choose_nonlinearity2(model_typ);
  for (int i = 0; i < S; i++) {
    mu_power_ncv.row(i) = model_fun(R_power_ncv.row(i), loudIn.doses);
  }

  returned_pval = pivotal_pvalue(R_power_ncv, &loudIn, mu_power_ncv);
  expected_pval = 0.001709650;
  expect_true(essentiallyEqual(expected_pval, returned_pval, 0.001));

  // Individual LogCV
  loudIn.dist = distribution::log_normal;
  loudIn.model = cont_model::exp_3;

  Eigen::MatrixXd R_exp3_log{
      {2.074359, 1.075104, 2.342920, 1.049676},
      {2.074831, 1.087564, 2.408487, 1.019421},
      {2.074691, 1.083938, 2.393829, 1.028519},
      {2.074695, 1.083803, 2.395533, 1.027121},
      {2.073364, 1.049481, 2.213176, 1.066825}
  };

  Eigen::MatrixXd mu_exp3_log(S, N);
  model_typ = getLoudModelType(loudIn.model, loudIn.dist, loudIn.datatype);
  model_fun = choose_nonlinearity2(model_typ);
  for (int i = 0; i < S; i++) {
    mu_exp3_log.row(i) = model_fun(R_exp3_log.row(i), loudIn.doses);
  }

  returned_pval = pivotal_pvalue(R_exp3_log, &loudIn, mu_exp3_log);
  expected_pval = 1.0;
  expect_true(essentiallyEqual(expected_pval, returned_pval, 0.001));

  // summary data
  loudIn.model = cont_model::power;
  Eigen::MatrixXd Y_sum{
      {10.61764, 0.8937421, 20},
      {11.54771, 1.0580638, 20},
      {12.20492, 1.3528275, 20},
      {14.73715, 1.0778448, 19},
      {15.85227, 0.8350199, 19}
  };

  Eigen::MatrixXd Dose_sum{{0}, {0.125}, {0.25}, {0.5}, {1.0}};

  Eigen::MatrixXd R_power_cv_sum{
      {11.71532, 15.99155, 1.593011, 0.5574789},
      {11.67326, 16.13868, 1.607036, 0.5692751},
      {11.66204, 16.17652, 1.610603, 0.5723202},
      {11.64077, 16.24593, 1.617101, 0.5779315},
      {11.63452, 16.26341, 1.618728, 0.5793729}
  };

  expected_pval = 1.0;

  loudIn.Y = Y_sum;
  loudIn.doses = Dose_sum;
  loudIn.datatype = loud_datatype::l_summary;

  // Summary CV
  loudIn.dist = distribution::normal;

  N = loudIn.Y.rows();
  Eigen::MatrixXd mu_power_cv_sum(S, N);
  model_typ = getLoudModelType(loudIn.model, loudIn.dist, loudIn.datatype);
  model_fun = choose_nonlinearity2(model_typ);
  for (int i = 0; i < S; i++) {
    mu_power_cv_sum.row(i) = model_fun(R_power_cv_sum.row(i), loudIn.doses);
  }

  returned_pval = pivotal_pvalue(R_power_cv_sum, &loudIn, mu_power_cv_sum);
  expect_true(essentiallyEqual(expected_pval, returned_pval, 0.001));

  // Summary NCV
  loudIn.dist = distribution::normal_ncv;
  Eigen::MatrixXd R_power_ncv_sum{
      {10.53703, 14.76691, 2.395046, 0.2182275, 1.257413},
      {10.53941, 14.76425, 2.395757, 0.2194998, 1.259096},
      {10.53811, 14.76550, 2.395344, 0.2188017, 1.258193},
      {10.54633, 14.75680, 2.397978, 0.2232843, 1.264136},
      {10.52846, 14.77626, 2.392139, 0.2132897, 1.250982}
  };

  Eigen::MatrixXd mu_power_ncv_sum(S, N);
  model_typ = getLoudModelType(loudIn.model, loudIn.dist, loudIn.datatype);
  model_fun = choose_nonlinearity2(model_typ);
  for (int i = 0; i < S; i++) {
    mu_power_ncv_sum.row(i) = model_fun(R_power_ncv_sum.row(i), loudIn.doses);
  }

  expected_pval = 0.0009404825;
  returned_pval = pivotal_pvalue(R_power_ncv_sum, &loudIn, mu_power_ncv_sum);
  expect_true(essentiallyEqual(expected_pval, returned_pval, 0.001));

  // Summary LogCV
  loudIn.dist = distribution::log_normal;
  loudIn.model = cont_model::exp_3;
  Eigen::MatrixXd R_exp3_log_sum{
      {3.378570, 1.876914, 0.9280766, 0.2711137},
      {3.376712, 1.891313, 1.0694113, 0.2856840},
      {3.376792, 1.890697, 1.0538501, 0.2841253},
      {3.377071, 1.888535, 1.0393825, 0.2826644},
      {3.377061, 1.888611, 1.0338722, 0.2821067}
  };

  Eigen::MatrixXd mu_exp3_log_sum(S, N);
  model_typ = getLoudModelType(loudIn.model, loudIn.dist, loudIn.datatype);
  model_fun = choose_nonlinearity2(model_typ);
  for (int i = 0; i < S; i++) {
    mu_exp3_log_sum.row(i) = model_fun(R_exp3_log_sum.row(i), loudIn.doses);
  }

  expected_pval = 1.0;
  returned_pval = pivotal_pvalue(R_exp3_log_sum, &loudIn, mu_exp3_log_sum);
  expect_true(essentiallyEqual(expected_pval, returned_pval, 0.001));
}

void rg_dg_test() {
  // rg
  int iter = 5000;
  //  Eigen::VectorXd log_mu {{11.0558576, 15.5372800, 0.3687901, -0.8283089}};
  //  Eigen::MatrixXd log_cov {
  //	  { 7.978730e-12, -3.390675e-08, -2.171972e-09, -2.115609e-07},
  //	  {-3.390675e-08,  1.440916e-04,  9.230106e-06,  8.990582e-04},
  //	  {-2.171972e-09,  9.230106e-06,  5.912550e-07,  5.759142e-05},
  //	  {-2.115609e-07,  8.990582e-04,  5.759142e-05,  5.612043e-03}
  //
  //  };
  Eigen::VectorXd log_mu{{10.1447911, 17.5885148, 0.2291463, -0.5226341}};
  Eigen::MatrixXd log_cov{
      {2.258106e-11, -7.429725e-11, 1.092654e-12, -1.905617e-09},
      {-7.429725e-11, 2.444563e-10, -3.595100e-12, 6.269948e-09},
      {1.092654e-12, -3.595100e-12, 5.287141e-14, -9.220911e-11},
      {-1.905617e-09, 6.269948e-09, -9.220911e-11, 1.608151e-07}
  };

  Eigen::MatrixXd g_estimate(iter, log_mu.size());
  std::vector<bool> isNegative(log_mu.size());
  isNegative[0] = true;
  isNegative[1] = true;
  for (int i = 2; i < isNegative.size(); i++) {
    isNegative[i] = false;
  }

  rg(iter, log_mu, log_cov, isNegative, g_estimate);
  Eigen::VectorXd means = g_estimate.colwise().mean().transpose();

  Eigen::VectorXd expectedMeans{{10.1447911, 17.5885146, 1.2575260, 0.5929534}};
  expect_true(expectedMeans.isApprox(means, 1.5e-6));

  // simple dmvnorm test
  Eigen::MatrixXd x{{0, 0}, {0, 0}, {0, 0}};
  Eigen::MatrixXd sigma{{4, 2}, {2, 3}};
  Eigen::VectorXd mu{{1, 2}};
  Eigen::VectorXd ret(x.rows());
  dmvnorm(x, mu, sigma, ret, true);

  Eigen::VectorXd expectedRes{{-3.5651, -3.5651, -3.5651}};
  expect_true(expectedRes.isApprox(ret, 1.5e-6));

  // dg test

  double expectedMean = 78.4598;
  Eigen::VectorXd ret2(g_estimate.rows());
  dg(g_estimate, log_mu, log_cov, isNegative, ret2);
  double mean = ret2.mean();
  expect_true(essentiallyEqual(expectedMean, mean, 0.001));
}

void bridge_sample_test() {
  // individual data
  Eigen::MatrixXd Y{{11.7692}, {11.5889}, {10.9485}, {12.0335}, {10.5843}, {10.9571}, {12.1292},
                    {12.5896}, {12.9235}, {11.1651}, {9.2031},  {11.7029}, {11.9342}, {11.8759},
                    {10.6886}, {12.0227}, {11.277},  {10.811},  {14.2668}, {10.4833}, {11.4515},
                    {12.5341}, {12.4558}, {13.7429}, {13.0073}, {12.5508}, {10.4007}, {12.8715},
                    {13.4137}, {9.33564}, {10.9298}, {10.6283}, {12.9618}, {12.2369}, {14.8123},
                    {11.0151}, {13.5539}, {10.7161}, {12.7046}, {12.7758}, {14.4463}, {14.8051},
                    {17.1731}, {14.4799}, {14.8711}, {13.9808}, {12.8826}, {15.3315}, {14.9358},
                    {14.0363}, {14.8257}, {14.0969}, {15.5807}, {15.202},  {14.1877}, {14.5575},
                    {13.1242}, {17.0857}, {14.4028}};

  Eigen::MatrixXd D{{0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125},
                    {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125}, {0.125},
                    {0.125}, {0.125}, {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},
                    {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},  {0.25},
                    {0.25},  {0.25},  {0.25},  {0.25},  {0.5},   {0.5},   {0.5},   {0.5},   {0.5},
                    {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},   {0.5},
                    {0.5},   {0.5},   {0.5},   {0.5},   {0.5}};

  // Individual CV
  Eigen::MatrixXd R{
      {10.14480, 17.58850, 1.257526, 0.5927195},
      {10.14479, 17.58852, 1.257526, 0.5929868},
      {10.14479, 17.58851, 1.257526, 0.5928710},
      {10.14479, 17.58851, 1.257526, 0.5928589},
      {10.14478, 17.58854, 1.257525, 0.5933468}
  };

  std::vector<double> mean = {10.6176, 11.5477, 12.2049, 14.7371, 15.8523};
  std::vector<int> N_obs = {20, 20, 20, 19, 19};
  std::vector<double> var = {0.798775, 1.1195, 1.83014, 1.16175, 0.697258};
  double ssq01 = 0.749389;
  int iter = 5;
  int burnin = BMDS_MISSING;  // Not used for this test, but needed for constructor function below
  int bmr = 1;
  int BMD_type = 2;
  bool isIncreasing = true;
  int weightOption = 3;  // WAIC & INT_FACTOR

  struct fitInput cvInput = createFitInput(
      D, Y, mean[0], mean[mean.size() - 1], N_obs[0], N_obs[N_obs.size() - 1], var[0],
      var[var.size() - 1], N_obs[0] + N_obs[N_obs.size() - 1], ssq01, iter, burnin, bmr,
      distribution::normal, loud_datatype::l_individual, BMD_type, isIncreasing, weightOption,
      BMDS_MISSING
  );

  cvInput.model = cont_model::power;
  Eigen::MatrixXd priorr{
      {5, 19.0000000, 10.61764, 0.2050385, 1},
      {5, 18.0000000, 15.85227, 0.1968161, 1},
      {2, 0.4700036, 0.42100, BMDS_MISSING, 1},
      {4, 19.0000000, 14.61307, BMDS_MISSING, 1}
  };

  struct fitResult cvPowerOut;
  std::vector<bool> isNegative(priorr.rows());
  // Eigen::VectorXd init(priorr.rows());
  // Eigen::MatrixXd diag = Eigen::VectorXd::Constant(priorr.rows(), 1.0).asDiagonal();
  isNegative[0] = true;
  isNegative[1] = true;
  for (int i = 2; i < isNegative.size(); i++) {
    isNegative[i] = false;
  }

  int model_typ = getLoudModelType(cvInput.model, cvInput.dist, cvInput.datatype);
  ptr2 model_fun = choose_nonlinearity2(model_typ);
  int S = R.rows();
  int N = cvInput.Y.rows();
  Eigen::MatrixXd mu(S, N);
  for (int i = 0; i < S; i++) {
    mu.row(i) = model_fun(R.row(i), cvInput.doses);
  }

  bridge_sample(R, mu, &cvInput, &cvPowerOut, priorr, isNegative);
  double expected_waic = -112.5453;
  expect_true(essentiallyEqual(expected_waic, cvPowerOut.waic, 0.001));

  double expected_intFactor = -202.7837;
  std::cout << "expected_intFactor:" << expected_intFactor
            << ", result int_factor:" << cvPowerOut.int_factor << std::endl;
  expect_true(essentiallyEqual(expected_intFactor, cvPowerOut.int_factor, 0.001));
  std::cout << "expected int_factor:" << expected_intFactor << std::endl;
  std::cout << "actual int_factor:" << cvPowerOut.int_factor << std::endl;
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
    // essentiallyEqual(expProbs1[i], probs[i], 1.5e-6);
    std::cout << "comparing " << expProbs1[i] << " to " << probs[i] << std::endl;
    expect_true(essentiallyEqual(expProbs1[i], probs[i], 0.001));
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
    expect_true(essentiallyEqual(expProbs1[i], probs[i], 0.001));
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

  expect_true(essentiallyEqual(lk, exp_lk, 0.001));
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
  expect_true(essentiallyEqual(ineqVal, expVal, 0.001));
  std::cout << "expected: " << expVal << std::endl;
  std::cout << "returned: " << ineqVal << std::endl;

  for (int i = 0; i < grad.size(); i++) {
    expect_true(essentiallyEqual(grad[i], expGrad[i], 0.001));
    std::cout << "expected: " << expGrad[i] << std::endl;
    std::cout << "returned: " << grad[i] << std::endl;
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

  expect_true(essentiallyEqual(eqVal, expVal, 0.001));

  for (int i = 0; i < grad.size(); i++) {
    expect_true(essentiallyEqual(grad[i], expGrad[i], 0.001));
  }
}

void cont_AIC_penalty_test() {
  struct continuous_analysis anal;
  anal.parms = 4;
  double prior[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1, 0.2, 1, -100, -100, 1, -18, 100, 100, 18, 18};
  anal.prior = prior;
  anal.model = power;

  struct BMDS_results bmdsRes;
  struct continuous_model_result res;
  double parm[4] = {5.05487, -0.0260177, 1, 0.698018};
  res.nparms = 4;
  res.parms = parm;
  res.max = 175.027;

  bool countAllParmsOnBoundary;
  countAllParmsOnBoundary = true;
  calcContAIC(&anal, &res, &bmdsRes, countAllParmsOnBoundary);
  int numBounded = 0;
  for (int i = 0; i < bmdsRes.bounded.size(); i++) {
    if (bmdsRes.bounded[i]) numBounded++;
  }

  double AIC_penalized = bmdsRes.AIC;

  countAllParmsOnBoundary = false;
  calcContAIC(&anal, &res, &bmdsRes, countAllParmsOnBoundary);
  double AIC_unpenalized = bmdsRes.AIC;

  expect_true(essentiallyEqual(AIC_unpenalized + 2 * numBounded, AIC_penalized, 0.001));
}

void dicho_AIC_penalty_test() {
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

  bool countAllParmsOnBoundary;
  countAllParmsOnBoundary = true;
  calcDichoAIC(&anal, &res, &bmdsRes, countAllParmsOnBoundary);
  int numBounded = 0;
  for (int i = 0; i < bmdsRes.bounded.size(); i++) {
    if (bmdsRes.bounded[i]) numBounded++;
  }

  double AIC_penalized = bmdsRes.AIC;

  countAllParmsOnBoundary = false;
  calcDichoAIC(&anal, &res, &bmdsRes, countAllParmsOnBoundary);
  double AIC_unpenalized = bmdsRes.AIC;

  expect_true(essentiallyEqual(AIC_unpenalized + 2 * numBounded, AIC_penalized, 0.001));
}

void nested_AIC_penalty_test() {
  double fitted_LL = -269.735;
  double fitted_df_pen = 30;
  double fitted_df_unpen = 35;
  double red_df = 38;
  int numBounded = fitted_df_unpen - fitted_df_pen;

  double AIC_penalized = calcNestedAIC(fitted_LL, fitted_df_pen, red_df);
  double AIC_unpenalized = calcNestedAIC(fitted_LL, fitted_df_unpen, red_df);

  expect_true(essentiallyEqual(AIC_unpenalized + 2 * numBounded, AIC_penalized, 0.001));
}

void additional_cont_calcs_test() {}

void additional_dicho_calcs_test() {
  struct dichotomous_analysis anal;
  struct dichotomous_model_result res;
  struct dichotomous_GOF gof;
  struct BMDS_results bmdsRes;
  struct dicho_AOD bmdsAOD;
  bool countAllParmsOnBoundary = false;
  bool isLoud = false;

  anal.model = 8;
  anal.n = 5;
  double Y[] = {0, 5, 30, 65, 90};
  double doses[] = {0, 50, 100, 150, 200};
  double n_group[] = {100, 100, 100, 100, 100};
  anal.prior_cols = 5;
  double prior[] = {0, 0, 0, 0, 0, 0, -18, 0, 18, 100};
  anal.BMD_type = 1;
  anal.BMR = 0.1;
  anal.alpha = 0.05;
  anal.degree = 2;
  anal.samples = 0;
  anal.burnin = 0;
  anal.parms = 2;
  anal.Y = Y;
  anal.doses = doses;
  anal.n_group = n_group;
  anal.prior = prior;

  res.model = 8;
  res.nparms = 2;
  double parms[] = {-18, 0.00595948};
  res.parms = parms;
  double cov[] = {3.03852e-09, -1.89595e-18, -1.89595e-18, 2.02955e-07};
  res.cov = cov;
  res.max = 210.797;
  res.dist_numE = 200;
  res.model_df = 2;
  res.total_df = 0;
  res.bmd = 17.6795;
  res.gof_p_value = 0;
  res.gof_chi_sqr_statistic = 0;
  double bmd_dist[] = {
      0,       0,       0,       15.0551, 15.1706, 15.2835, 15.3813,  15.4611,  15.5282,  15.5877,
      15.6446, 15.7011, 15.7525, 15.7992, 15.8423, 15.8829, 15.9219,  15.9603,  15.9988,  16.0355,
      16.0703, 16.1034, 16.1351, 16.1655, 16.195,  16.2237, 16.252,   16.28,    16.3081,  16.3357,
      16.3625, 16.3885, 16.4137, 16.4383, 16.4623, 16.4859, 16.509,   16.5317,  16.5541,  16.5764,
      16.5985, 16.6205, 16.6425, 16.6643, 16.6857, 16.7067, 16.7274,  16.7478,  16.7679,  16.7878,
      16.8074, 16.8268, 16.846,  16.8651, 16.884,  16.9029, 16.9217,  16.9404,  16.9591,  16.9778,
      16.9964, 17.0148, 17.0331, 17.0512, 17.0692, 17.087,  17.1047,  17.1223,  17.1397,  17.1571,
      17.1744, 17.1916, 17.2088, 17.226,  17.2431, 17.2602, 17.2772,  17.2943,  17.3114,  17.3286,
      17.3456, 17.3626, 17.3795, 17.3963, 17.4131, 17.4298, 17.4465,  17.4631,  17.4797,  17.4963,
      17.5128, 17.5294, 17.5459, 17.5625, 17.5791, 17.5957, 17.6123,  17.629,   17.6458,  17.6626,
      17.6795, 17.6964, 17.7133, 17.7301, 17.7469, 17.7638, 17.7806,  17.7974,  17.8143,  17.8312,
      17.8481, 17.8651, 17.8821, 17.8992, 17.9164, 17.9337, 17.9511,  17.9685,  17.9861,  18.0039,
      18.0217, 18.0397, 18.0577, 18.0758, 18.0938, 18.1119, 18.1301,  18.1483,  18.1666,  18.185,
      18.2035, 18.2222, 18.241,  18.2599, 18.2791, 18.2984, 18.3179,  18.3377,  18.3577,  18.3779,
      18.3984, 18.419,  18.4396, 18.4603, 18.4811, 18.502,  18.5231,  18.5444,  18.5659,  18.5877,
      18.6098, 18.6323, 18.6551, 18.6783, 18.702,  18.7262, 18.7508,  18.7759,  18.8011,  18.8264,
      18.8519, 18.8776, 18.9037, 18.9303, 18.9575, 18.9853, 19.0138,  19.0431,  19.0734,  19.1046,
      19.137,  19.1698, 19.2028, 19.2361, 19.27,   19.3049, 19.341,   19.3785,  19.4179,  19.4594,
      19.5031, 19.5489, 19.5947, 19.6414, 19.69,   19.7416, 19.7972,  19.8579,  19.9244,  19.9922,
      20.0622, 20.1391, 20.2273, 20.3313, 20.4541, 20.587,  INFINITY, INFINITY, INFINITY, INFINITY,
      0,       0.005,   0.01,    0.015,   0.02,    0.025,   0.03,     0.035,    0.04,     0.045,
      0.05,    0.055,   0.06,    0.065,   0.07,    0.075,   0.08,     0.085,    0.09,     0.095,
      0.1,     0.105,   0.11,    0.115,   0.12,    0.125,   0.13,     0.135,    0.14,     0.145,
      0.15,    0.155,   0.16,    0.165,   0.17,    0.175,   0.18,     0.185,    0.19,     0.195,
      0.2,     0.205,   0.21,    0.215,   0.22,    0.225,   0.23,     0.235,    0.24,     0.245,
      0.25,    0.255,   0.26,    0.265,   0.27,    0.275,   0.28,     0.285,    0.29,     0.295,
      0.3,     0.305,   0.31,    0.315,   0.32,    0.325,   0.33,     0.335,    0.34,     0.345,
      0.35,    0.355,   0.36,    0.365,   0.37,    0.375,   0.38,     0.385,    0.39,     0.395,
      0.4,     0.405,   0.41,    0.415,   0.42,    0.425,   0.43,     0.435,    0.44,     0.445,
      0.45,    0.455,   0.46,    0.465,   0.47,    0.475,   0.48,     0.485,    0.49,     0.495,
      0.5,     0.505,   0.51,    0.515,   0.52,    0.525,   0.53,     0.535,    0.54,     0.545,
      0.55,    0.555,   0.56,    0.565,   0.57,    0.575,   0.58,     0.585,    0.59,     0.595,
      0.6,     0.605,   0.61,    0.615,   0.62,    0.625,   0.63,     0.635,    0.64,     0.645,
      0.65,    0.655,   0.66,    0.665,   0.67,    0.675,   0.68,     0.685,    0.69,     0.695,
      0.7,     0.705,   0.71,    0.715,   0.72,    0.725,   0.73,     0.735,    0.74,     0.745,
      0.75,    0.755,   0.76,    0.765,   0.77,    0.775,   0.78,     0.785,    0.79,     0.795,
      0.8,     0.805,   0.81,    0.815,   0.82,    0.825,   0.83,     0.835,    0.84,     0.845,
      0.85,    0.855,   0.86,    0.865,   0.87,    0.875,   0.88,     0.885,    0.89,     0.895,
      0.9,     0.905,   0.91,    0.915,   0.92,    0.925,   0.93,     0.935,    0.94,     0.945,
      0.95,    0.955,   0.96,    0.965,   0.97,    0.975,   0.98,     0.985,    0.99,     0.995,
  };
  res.bmd_dist = bmd_dist;

  bmdsRes.bounded.push_back(false);
  bmdsRes.bounded.push_back(false);

  additional_dicho_calcs(&anal, &res, &gof, &bmdsRes, &bmdsAOD, &countAllParmsOnBoundary, &isLoud);

  double expBMD = 17.6795;
  double expBMDL = 15.6446;
  double expBMDU = 20.0622;
  expect_true(essentiallyEqual(bmdsRes.BMD, expBMD, 1e-4));
  expect_true(essentiallyEqual(bmdsRes.BMDL, expBMDL, 1e-4));
  expect_true(essentiallyEqual(bmdsRes.BMDU, expBMDU, 1e-4));

  std::vector<double> gofExpected = {1.523e-06, 25.7679, 44.896, 59.0952, 69.6355};
  std::vector<double> gofRes = {-0.0012341, -4.74851, -2.99485, 1.201, 4.42869};
  for (int i = 0; i < gofExpected.size(); i++) {
    expect_true(essentiallyEqual(gof.expected[i], gofExpected[i], 1e-4));
    expect_true(essentiallyEqual(gof.residual[i], gofRes[i], 1e-4));
  }

  double expAIC = 423.594;
  double expBIC = -195.124;
  double expChisq = 52.5732;

  expect_true(essentiallyEqual(bmdsRes.AIC, expAIC, 1e-3));
  expect_true(essentiallyEqual(bmdsRes.BIC_equiv, expBIC, 1e-3));
  expect_true(essentiallyEqual(bmdsRes.chisq, expChisq, 1e-3));

  double expFullLL = -178.191;
  int expNFull = 5;
  double expRedLL = -332.032;
  int expNRed = 1;
  double expFitLL = -210.797;
  int expNFit = 1;
  double expDevFit = 65.212;
  double expDevRed = 307.682;
  int expDFFit = 4;
  int expDFRed = 4;
  double expPVFit = 2.32148e-13;
  double expPVRed = 0;

  expect_true(essentiallyEqual(bmdsAOD.fullLL, expFullLL, 1e-3));
  expect_true(essentiallyEqual(bmdsAOD.nFull, expNFull, 1));
  expect_true(essentiallyEqual(bmdsAOD.redLL, expRedLL, 1e-3));
  expect_true(essentiallyEqual(bmdsAOD.nRed, expNRed, 1));
  expect_true(essentiallyEqual(bmdsAOD.fittedLL, expFitLL, 1e-3));
  expect_true(essentiallyEqual(bmdsAOD.nFit, expNFit, 1));
  expect_true(essentiallyEqual(bmdsAOD.devFit, expDevFit, 1e-3));
  expect_true(essentiallyEqual(bmdsAOD.devRed, expDevRed, 1e-3));
  expect_true(essentiallyEqual(bmdsAOD.dfFit, expDFFit, 1));
  expect_true(essentiallyEqual(bmdsAOD.dfRed, expDFRed, 1));
  expect_true(essentiallyEqual(bmdsAOD.pvFit, expPVFit, 1e-3));
  expect_true(essentiallyEqual(bmdsAOD.pvRed, expPVRed, 1e-3));
}
