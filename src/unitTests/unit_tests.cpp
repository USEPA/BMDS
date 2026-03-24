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
  nested_AIC_penalty_test();
  run_loud_model_fit_test();

  return 0;
}

void objfunc_test() {
  std::vector<double> x{1.5, 2.0, 3.2};
  std::vector<double> tmp;
  // assert(objfunc_bmdl(x, tmp, NULL)==1.5);
  expect_true(objfunc_bmdl(x, tmp, NULL) == 1.5);
}

// these tests expect alpha to be returned not ln(alpha) as in BMDS
void run_loud_model_fit_test() {
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
      {10.39001, 5.213014, 1.515874, 1.6183734, -2.8310805},
      {10.39067, 5.213711, 1.515864, 1.7689691, -3.1755334},
      {10.39022, 5.213237, 1.515871, 1.6539914, -2.9111840},
      {10.39023, 5.213244, 1.515872, 1.6917375, -2.9998844},
      {10.37943, 5.201839, 1.516049, 0.2849009, 0.1829478}
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
      {10.74145, 0.5789798, 1.487489, -4.158105, 11.03095},
      {10.74138, 0.5790119, 1.487592, -4.166388, 11.05439},
      {10.74145, 0.5789768, 1.487479, -4.157356, 11.02883},
      {10.74150, 0.5789575, 1.487417, -4.152342, 11.01464},
      {10.74081, 0.5792735, 1.488431, -4.235545, 11.25011}
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
      {7.959443, 0.9996822, 2.342920, 0.9526748},
      {7.963203, 0.9946935, 2.408487, 0.9809492},
      {7.962085, 0.9961268, 2.393829, 0.9722722},
      {7.962117, 0.9961875, 2.395533, 0.9735949},
      {7.951526, 1.0107213, 2.213176, 0.9373611}
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
      {10.81164, 1.973901, 1.559089, 2.442058, 0.3253742, -0.5308331},
      {10.83474, 1.969843, 1.557273, 2.433451, 0.2547470, -0.3016726},
      {10.82757, 1.971251, 1.557869, 2.436323, 0.2916227, -0.4211056},
      {10.82954, 1.971680, 1.557864, 2.436824, 0.2845722, -0.3982736},
      {10.78528, 1.979022, 1.561160, 2.452975, 0.4082274, -0.8039603}
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
      {5.748410, 2.564370, 3.456401, 1.0683093, 1.311286},
      {5.759851, 2.565635, 3.441258, 1.0621063, 1.313785},
      {5.790654, 2.568815, 3.401958, 1.0435225, 1.321313},
      {5.821355, 2.572834, 3.363387, 1.0251954, 1.328821},
      {5.924718, 2.585508, 3.241374, 0.9598605, 1.356212}
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
      {12.09806, 55.54020, 2.428470, 2.827084, 3.133312, -7.098091},
      {12.09780, 55.73431, 2.434746, 2.823307, 3.257137, -7.456318},
      {12.09793, 55.63348, 2.431686, 2.825010, 3.162207, -7.181550},
      {12.09854, 55.11642, 2.416787, 2.832798, 2.783641, -6.090638},
      {12.09880, 54.88379, 2.410370, 2.835964, 2.630579, -5.651274}
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
      {26.47598, 4.119359, -1.176881, 0.9196827, 0.4113099},
      {26.47217, 4.118264, -1.179743, 0.9215378, 0.4119523},
      {26.47202, 4.118268, -1.179761, 0.9215637, 0.4116585},
      {26.47263, 4.118428, -1.179314, 0.9212743, 0.4116270},
      {26.48417, 4.121061, -1.172655, 0.9165915, 0.4097042}
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
      {8.765806, 0.10712713, 3.369291, 4.113787, 0.3856143},
      {8.764970, 0.08819036, 3.328424, 4.104144, 0.3848478},
      {8.763489, 0.06188680, 3.274165, 4.091029, 0.3837737},
      {8.764894, 0.09084669, 3.334945, 4.105468, 0.3848634},
      {8.763918, 0.07261393, 3.296879, 4.096308, 0.3841136}
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
      {10.76265, 1.329088, 3.166795, 0.1593405, 12.93630},
      {10.76173, 1.329261, 3.168671, 0.1598616, 12.93521},
      {10.76033, 1.329521, 3.171500, 0.1609399, 12.93295},
      {10.75474, 1.330652, 3.183645, 0.1634249, 12.92774},
      {10.75480, 1.330088, 3.178540, 0.1626546, 12.92935}
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
  expect_true(log_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << log_efsaLogCV_reldev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

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
      {30.33905, 1.953997, 0.1791335, 0.5720145, 0.3390048},
      {30.33906, 1.953983, 0.1791344, 0.5719660, 0.3390056},
      {30.33910, 1.953895, 0.1791402, 0.5716567, 0.3390111},
      {30.33909, 1.953915, 0.1791389, 0.5717267, 0.3390099},
      {30.33910, 1.953901, 0.1791398, 0.5716806, 0.3390107}
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
      {13.00834, 1.294782, 1.173444, 2.019590, 3.057626},
      {13.00392, 1.299844, 1.172860, 2.021746, 3.048809},
      {13.05378, 1.244163, 1.179499, 1.999033, 3.150114},
      {13.04187, 1.257789, 1.177901, 2.004633, 3.124734},
      {12.99180, 1.313767, 1.171260, 2.027419, 3.024512}
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
  expect_true(lms_efsaLogCV_reldev_BMD.isApprox(fitOut.BMD, 1.5e-6));
  std::cout << "expected BMDs:" << std::endl << lms_efsaLogCV_reldev_BMD << std::endl;
  std::cout << "actual BMDs:" << std::endl << fitOut.BMD << std::endl;

  //  std::cout<<"expected parms:"<<std::endl<<exp5LogCVParms<<std::endl;
  //  std::cout<<"actual parms:"<<std::endl<<fitOut.parms<<std::endl;

  //  std::cout<<"expected BMDs:"<<std::endl<<exp5LogCV_stdev_BMD<<std::endl;
  //  std::cout<<"actual BMDs:"<<std::endl<<fitOut.BMD<<std::endl;
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

  double AIC_penalized = bmdsRes.AIC;

  countAllParmsOnBoundary = false;
  calcContAIC(&anal, &res, &bmdsRes, countAllParmsOnBoundary);
  double AIC_unpenalized = bmdsRes.AIC;

  essentiallyEqual(AIC_unpenalized - 2, AIC_penalized, 1e-6);
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

  double AIC_penalized = bmdsRes.AIC;

  countAllParmsOnBoundary = false;
  calcDichoAIC(&anal, &res, &bmdsRes, countAllParmsOnBoundary);
  double AIC_unpenalized = bmdsRes.AIC;

  essentiallyEqual(AIC_unpenalized - 2, AIC_penalized, 1e-6);
}

void nested_AIC_penalty_test() {
  double fitted_LL = -269.735;
  double fitted_df_pen = 30;
  double fitted_df_unpen = 35;
  double red_df = 38;
  int numBounded = 5;

  double AIC_penalized = calcNestedAIC(fitted_LL, fitted_df_pen, red_df);
  double AIC_unpenalized = calcNestedAIC(fitted_LL, fitted_df_unpen, red_df);

  essentiallyEqual(AIC_unpenalized - 2 * numBounded, AIC_penalized, 1e-6);
}
