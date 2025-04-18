#include <vector>
/////////////////////////
// Initial values & Priors//
///////////////////////////
//
// Dichotomous Models
//
//////////BAYESIAN////////////////
std::vector<double> prBayesianDHill{1,   1,   1,   2,   -1,  0, -3, 0.693147, 2,  3,
                                    3.3, 0.5, -40, -40, -40, 0, 40, 40,       40, 40};
std::vector<double> prBayesianGamma{1, 2,   2,   0, 0.693147, 0,  2,    0.424264,
                                    1, -18, 0.2, 0, 18,       20, 10000};
std::vector<double> prBayesianLogistic{1, 2, 0, 0, 2, 2, -20, 0, 20, 40};
std::vector<double> prBayesianLogLogistic{1,   1,   2,   0, 0,  0.693147, 2, 1,
                                          0.5, -20, -40, 0, 20, 40,       20};
std::vector<double> prBayesianLogProbit{1,   1,   2,   0, 0,  0.693147, 2, 1,
                                        0.5, -20, -40, 0, 20, 40,       20};
std::vector<double> prBayesianMulti1{1, 2, 0, 0, 2, 1, -20, 0, 20, 1e6};  // degree 1
std::vector<double> prBayesianMulti2{1, 2,   2, 0, 0,  0,   2,  1,
                                     1, -20, 0, 0, 20, 1e6, 1e6};  // degree 2
std::vector<double> prBayesianMulti3{1, 2, 2,   2, 0, 0, 0,  0,   2,   1,
                                     1, 1, -20, 0, 0, 0, 20, 1e6, 1e6, 1e6};  // degree 3
std::vector<double> prBayesianMulti4{1, 2, 2,   2, 2, 0, 0, 0,  0,   0,   2,   1,  1,
                                     1, 1, -20, 0, 0, 0, 0, 20, 1e6, 1e6, 1e6, 1e6};  // degree 4
std::vector<double> prBayesianMulti5{
    1, 2, 2, 2,   2, 2, 0, 0, 0, 0,  0,   0,   2,   1,   1,
    1, 1, 1, -20, 0, 0, 0, 0, 0, 20, 1e6, 1e6, 1e6, 1e6, 1e6
};  // degree 5
std::vector<double> prBayesianProbit{1, 2, 0, 0, 2, 1, -20, 0, 20, 40};
std::vector<double> prBayesianQLinear{1, 2, 0, 0, 2, 1, -20, 0, 20, 18};
std::vector<double> prBayesianWeibull{1,   2,   2, 0, 0.424264, 0,  2,    0.5,
                                      1.5, -20, 0, 0, 20,       40, 10000};
//////////UNRESTRICTED FREQ////////////////
std::vector<double> prUFreqLogistic{0, 0, 0, 0, 0, 0, -18, 0, 18, 100};
std::vector<double> prUFreqDHill{0, 0, 0,   0,   0,   0,    0,  0,  0,  0,
                                 0, 0, -18, -18, -18, 1e-8, 18, 18, 18, 18};
std::vector<double> prUFreqGamma{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 0.2, 0, 18, 18, 100};
std::vector<double> prUFreqLogLogistic{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, -18, 1e-4, 18, 18, 18};
std::vector<double> prUFreqLogProbit{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, -18, 1e-4, 18, 18, 18};
std::vector<double> prUFreqMulti1{0, 0, 0, 0, 0, 0, -18, -18, 18, 100};  // degree 1
std::vector<double> prUFreqMulti2{0, 0,   0,   0,   0,  0,   0,  0,
                                  0, -18, -18, -18, 18, 100, 1e4};  // degree 2
std::vector<double> prUFreqMulti3{0, 0, 0,   0,   0,   0,   0,  0,   0,   0,
                                  0, 0, -18, -18, -18, -18, 18, 100, 1e4, 1e4};  // degree 3
std::vector<double> prUFreqMulti4{
    0, 0, 0,   0,   0,   0,   0,   0,  0,   0,   0,   0,  0,
    0, 0, -18, -18, -18, -18, -18, 18, 100, 1e4, 1e4, 1e4
};  // degree 4 NOTWORKING
std::vector<double> prUFreqMulti5{0,   0,   0,  0,   0,   0,   0,   0,   0,   0,   0,
                                  0,   0,   0,  0,   0,   0,   0,   -18, -18, -18, -18,
                                  -18, -18, 18, 1e4, 1e4, 1e4, 1e4, 1e4};  // degree 5 NOTWORKING
std::vector<double> prUFreqProbit{0, 0, 0, 0, 0, 0, -18, 0, 18, 18};
std::vector<double> prUFreqQLinear{0, 0, 0, 0, 0, 0, -18, 0, 18, 100};
std::vector<double> prUFreqWeibull{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 1e-6, 1e-6, 18, 18, 100};
// RESTRICTED FREQ////////////////
std::vector<double> prRFreqDHill{0, 0, 0,   0,   0,   0, 0,  0,  0,  0,
                                 0, 0, -18, -18, -18, 1, 18, 18, 18, 18};
std::vector<double> prRFreqGamma{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 1, 0, 18, 18, 100};
std::vector<double> prRFreqLogLogistic{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, -18, 1, 18, 18, 18};
std::vector<double> prRFreqLogProbit{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, -18, 1, 18, 18, 18};
std::vector<double> prRFreqMulti1{0, 0, 0, 0, 0, 0, -18, 0, 18, 1e4};                   // degree 1
std::vector<double> prRFreqMulti2{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 0, 0, 18, 1e4, 1e4};  // degree 2
std::vector<double> prRFreqMulti3{0, 0, 0,   0, 0, 0, 0,  0,   0,   0,
                                  0, 0, -18, 0, 0, 0, 18, 1e4, 1e4, 1e4};  // degree 3
std::vector<double> prRFreqMulti4{0, 0, 0,   0, 0, 0, 0, 0,  0,   0,   0,   0,  0,
                                  0, 0, -18, 0, 0, 0, 0, 18, 1e4, 1e4, 1e4, 1e4};  // degree 4
std::vector<double> prRFreqMulti5{
    0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,
    0, 0, 0, -18, 0, 0, 0, 0, 0, 18, 1e4, 1e4, 1e4, 1e4, 1e4
};  // degree 5
std::vector<double> prRFreqWeibull{0, 0, 0, 0, 0, 0, 0, 0, 0, -18, 1, 1e-6, 18, 18, 100};
/////////////////////////////////////////////////////////////////////////////
// Continuous Models
/////////////////////////////////////////////////////////////////////////////
// BAYESIAN
std::vector<double> prBayesianHill{2,      1, 2, 2,   1, 0, 1,   -0.69315, 0.405465, 0,  1,  2, 1,
                                   0.2501, 1, 0, -18, 0, 0, -18, 18,       18,       18, 18, 18};
std::vector<double> prBayesianHillNCV{2, 1,   2,  2,  1,      0,  1, -0.69315, 0.405465,
                                      0, 1,   2,  1,  0.2501, 1,  0, -18,      0,
                                      0, -18, 18, 18, 18,     18, 18};
std::vector<double> prBayesianPower{2,   1, 2, 1,      0, 0,   0.405465, 0,   1,  1,
                                    0.5, 1, 0, -10000, 0, -18, 1e6,      1e4, 40, 18};
std::vector<double> prBayesianPowerNCV{2, 1,   2,   2,   1,      0,  0, 0.405465, 0,
                                       0, 1,   1,   0.5, 0.2501, 1,  0, -10000,   0,
                                       0, -18, 1e6, 1e4, 40,     18, 18};
// funl
std::vector<double> prBayesianPoly1{2, 1, 1, 0, 0, 0, 1, 2, 1, 0, -10000, -18, 1e6, 1e4, 18};
std::vector<double> prBayesianPoly2{2, 1, 1, 1,      0,      0,   0,   0,   1,   2,
                                    2, 1, 0, -10000, -10000, -18, 1e6, 1e4, 1e4, 18};  // poly 2
std::vector<double> prBayesianPoly3{
    2, 1, 1, 1,      1,      0,      0,   0,   0,   0,   1,   2, 2,
    2, 1, 0, -10000, -10000, -10000, -18, 1e6, 1e4, 1e4, 1e4, 18
};  // poly 3
std::vector<double> prBayesianPoly4{2,      1,   1,   1,   1,   1,   0,   0, 0,      0,      0,
                                    0,      1,   2,   2,   2,   1,   1,   0, -10000, -10000, -10000,
                                    -10000, -18, 1e6, 1e4, 1e4, 1e4, 1e4, 18};  // poly 4
std::vector<double> prBayesianPoly5{2,   1,   1,   1,   1,      1,      1,      0,      0,
                                    0,   0,   0,   0,   0,      1,      2,      2,      2,
                                    2,   1,   1,   0,   -10000, -10000, -10000, -10000, -10000,
                                    -18, 1e6, 1e4, 1e4, 1e4,    1e4,    1e4,    18};  // poly 5
std::vector<double> prBayesianPoly1NCV{2,      1, 2, 1,      0, 0,   0,   0,   1,  2,
                                       0.2501, 1, 0, -10000, 0, -18, 1e6, 1e4, 18, 18};
std::vector<double> prBayesianPoly2NCV{2, 1,   1,   2,   1,      0,  0, 0,      0,
                                       0, 1,   2,   2,   0.2501, 1,  0, -10000, -10000,
                                       0, -18, 1e6, 1e4, 1e4,    18, 18};  // poly 2
std::vector<double> prBayesianPoly3NCV{2, 1,   1,   1,   2,   1,      0,  0, 0,      0,      0,
                                       0, 1,   2,   2,   2,   0.2501, 1,  0, -10000, -10000, -10000,
                                       0, -18, 1e6, 1e4, 1e4, 1e4,    18, 18};  // poly 3
std::vector<double> prBayesianPoly4NCV{
    2, 1,      1, 1, 1,      2,      1,      0,      0, 0,   0,   0,   0,   0,   1,   2,  2, 2,
    2, 0.2501, 1, 0, -10000, -10000, -10000, -10000, 0, -18, 1e6, 1e4, 1e4, 1e4, 1e4, 18, 18
};  // poly 4
std::vector<double> prBayesianPoly5NCV{2,      1,      1,      1, 1,      1,   2,   1,      0,
                                       0,      0,      0,      0, 0,      0,   0,   1,      2,
                                       2,      2,      2,      2, 0.2501, 1,   0,   -10000, -10000,
                                       -10000, -10000, -10000, 0, -18,    1e6, 1e4, 1e4,    1e4,
                                       1e4,    1e4,    18,     18};  // poly 5
std::vector<double> prBayesianExp5{2,      2, 1, 2, 1,   0, 0,   0,   0,   0,  1,  1, 1,
                                   0.2501, 1, 0, 0, -20, 0, -18, 1e6, 100, 20, 18, 18};
std::vector<double> prBayesianExp5NCV{2,   2, 1, 2,   2,   1,      0,   0,  0,  0,
                                      0,   0, 1, 1,   1,   0.2501, 0.5, 1,  0,  0,
                                      -20, 0, 0, -18, 1e6, 100,    20,  18, 18, 18};
// UNRESTRICTED FREQ
std::vector<double> prUFreqHillNormal{
    0,   0, 0,    0,    0, 0,    0,   0,   1,   0, 1,  2, 1,
    1.2, 1, -1e2, -1e2, 0, 1e-8, -18, 1e2, 1e2, 5, 18, 18
};  // normal dist
std::vector<double> prUFreqHillNormalNCV{0,    0,    0,   0,   0,  0,  0,    0,    0,    0, 0,
                                         0,    1,    1,   1,   0,  1,  1,    -1e8, -1e8, 0, 1e-8,
                                         -1e3, -1e3, 1e8, 1e8, 30, 18, 1000, 1000};  // normal dist
std::vector<double> prUFreqHillLognormal{0,    0,    0,   0,   0,   0,   0,    0,    0,
                                         0,    1,    1,   1,   0,   1,   1e-8, -1e8, 0,
                                         1e-8, -1e3, 1e8, 1e8, 100, 100, 1000};  // normal dist
std::vector<double> prUFreqPower{0,   0, 0,    0,    0, 0,   1,   0,   0.1, 1,
                                 0.2, 1, -100, -100, 0, -18, 100, 100, 18,  18};
std::vector<double> prUFreqPowerNCV{0,     0,     0,   0,   0,   0,    0,    0,    0,
                                    0,     1,     1,   1,   1,   1,    1e-8, -1e8, 1e-8,
                                    -1000, -1000, 1e8, 1e8, 100, 1000, 1000};
// funl
// priors for auto detect adv dir
std::vector<double> prUFreqPoly1{0, 0,    0,    0,   0,   0,   0, 0,
                                 0, -1e6, -1e6, -18, 1e6, 1e6, 18};  // poly 1
std::vector<double> prUFreqPoly2{0, 0, 0,    0,    0,    0,   0,   0,   0,   0,
                                 0, 0, -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 18};  // poly 2
std::vector<double> prUFreqPoly3{
    0, 0, 0,    0,    0,    0,    0,   0,   0,   0,   0,   0, 0,
    0, 0, -1e6, -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 1e6, 18
};  // poly 3
std::vector<double> prUFreqPoly4{0,    0,    0,    0,   0,   0,   0,   0,   0,    0,
                                 0,    0,    0,    0,   0,   0,   0,   0,   -1e6, -1e6,
                                 -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 1e6, 1e6,  18};  // poly 4
std::vector<double> prUFreqPoly5{
    0, 0, 0, 0,    0,    0,    0,    0,    0,    0,   0,   0,   0,   0,   0,   0,   0, 0,
    0, 0, 0, -1e6, -1e6, -1e6, -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 18
};  // poly 5
std::vector<double> prUFreqPoly1NCV{0, 0, 0, 0,   0, 0,   0,    0,  0,  0,
                                    0, 0, 0, -18, 0, -18, 1000, 18, 18, 18};  // poly 1
std::vector<double> prUFreqPoly2NCV{0, 0, 0, 0,   0,    0, 0,   0,    0,  0,   0,  0, 0,
                                    0, 0, 0, -18, -1e6, 0, -18, 1000, 18, 1e6, 18, 18};  // poly 2
std::vector<double> prUFreqPoly3NCV{0,    0,    0, 0,   0,    0,  0,   0,   0,  0,
                                    0,    0,    0, 0,   0,    0,  0,   0,   0,  -18,
                                    -1e6, -1e6, 0, -18, 1000, 18, 1e6, 1e6, 18, 18};  // poly 3
std::vector<double> prUFreqPoly4NCV{0,    0,    0, 0,   0,    0,  0,   0,   0,   0,  0,   0,
                                    0,    0,    0, 0,   0,    0,  0,   0,   0,   0,  -18, -1e6,
                                    -1e6, -1e6, 0, -18, 1000, 18, 1e6, 1e6, 1e6, 18, 18};  // poly 4
std::vector<double> prUFreqPoly5NCV{0, 0,   0,    0,  0,   0,   0,    0,    0,    0,
                                    0, 0,   0,    0,  0,   0,   0,    0,    0,    0,
                                    0, 0,   0,    0,  0,   -18, -1e6, -1e6, -1e6, -1e6,
                                    0, -18, 1000, 18, 1e6, 1e6, 1e6,  1e6,  18,   18};  // poly 5

// RESTRICTED FREQ
std::vector<double> prRFreqExp5Normal{0,   0, 0, 0, 0,   0, 0,   0,   0,   0,  0.1, 1, 0.5,
                                      0.2, 1, 0, 0, -20, 1, -18, 100, 100, 20, 18,  18};
std::vector<double> prRFreqExp5NormalNCV{0,   0, 0,   0,   0,   0,   0,   0,  0,  0,
                                         0,   0, 0.1, 1,   0.5, 0.2, 0.5, 1,  0,  0,
                                         -20, 1, 0,   -18, 100, 100, 20,  18, 18, 18};
std::vector<double> prRFreqExp5Lognormal{0,   0, 0,     0, 0,   0, 0,   0,    1,   0,  0.1, 1, 1,
                                         0.2, 1, -1000, 0, -20, 1, -18, 1000, 100, 20, 18,  18};
std::vector<double> prRFreqHillNormal{
    0,   0, 0,    0,    0, 0, 0,   0,   0,   0, 1,  2, 1,
    1.2, 1, -100, -100, 0, 1, -18, 100, 100, 5, 18, 18
};  // normal dist
std::vector<double> prRFreqHillNormalNCV{0,    0,    0,   0,    0,  0,  0,    0,    0,     0, 0,
                                         0,    1,    1,   1,    0,  1,  1,    -1e8, -1000, 0, 1,
                                         -1e3, -1e3, 1e8, 1000, 30, 18, 1000, 1000};  // normal dist
std::vector<double> prRFreqHillLognormal{
    0, 0, 0,    0,    0, 0, 0,    0,   0,   0,   1,   1,   1,
    0, 1, 1e-8, -1e8, 0, 1, -1e3, 1e8, 1e8, 100, 100, 1000
};  // normal dist
std::vector<double> prRFreqPower{0,   0, 0,    0,    0, 0,   0,   0,   0.1, 1,
                                 0.2, 1, -100, -100, 1, -18, 100, 100, 18,  18};  // SEG FAULT
std::vector<double> prRFreqPowerNCV{0,   0, 0,    0,      0, 0,   1,   0,     1,  1,
                                    0.2, 1, -100, -10000, 0, -18, 100, 10000, 18, 18};  // SEG FAULT
// funl
// priors for auto detect adv dir
std::vector<double> prRFreqPoly1{0, 0,     0,   0,   0,    0,  5, 5,
                                 1, -1000, -18, -18, 1000, 18, 18};  // poly 1
std::vector<double> prRFreqPoly2{0, 0, 0,    0,    0,    0,   0,   0,   5,   5,
                                 5, 1, -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 18};  // poly 2
std::vector<double> prRFreqPoly3{
    0, 0, 0,    0,    0,    0,    0,   0,   0,   0,   5,   5, 5,
    5, 1, -1e6, -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 1e6, 18
};  // poly 3
std::vector<double> prRFreqPoly4{0,    0,    0,    0,   0,   0,   0,   0,   0,    0,
                                 0,    0,    5,    5,   5,   5,   5,   1,   -1e6, -1e6,
                                 -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 1e6, 1e6,  18};  // poly 4
std::vector<double> prRFreqPoly5{
    0, 0, 0, 0,    0,    0,    0,    0,    0,    0,   0,   0,   0,   0,   5,   5,   5, 5,
    5, 5, 1, -1e6, -1e6, -1e6, -1e6, -1e6, -1e6, -18, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 18
};  // poly 5
std::vector<double> prRFreqPoly1NCV{0, 0, 0,    0,    0,    0,    0,   0,   1,    1,
                                    1, 1, -1e8, -1e8, 1000, -1e8, 1e8, 1e8, 1000, 1e8};  // poly 1
std::vector<double> prRFreqPoly2NCV{0, 0, 0, 0,   0,    0, 0,   0,    0,  0,   5,  5, 5,
                                    1, 1, 0, -18, -1e6, 0, -18, 1000, 18, 1e6, 18, 18};  // poly 2
std::vector<double> prRFreqPoly3NCV{0,    0,     0,     0,      0,      0,      0, 0,
                                    0,    0,     0,     0,      5,      5,      5, 5,
                                    1,    1,     -1000, -10000, -10000, -10000, 0, -18,
                                    1000, 10000, 10000, 10000,  100,    18};  // poly 3
std::vector<double> prRFreqPoly4NCV{0,    0,   0,   0,    0,    0,    0,    0,    0,
                                    0,    0,   0,   0,    0,    1,    1,    1,    1,
                                    1,    1,   1,   -1e8, -1e8, -1e8, -1e8, -1e8, -1000,
                                    -1e8, 1e8, 1e8, 1e8,  1e8,  1e8,  1000, 1e8};  // poly 4
std::vector<double> prRFreqPoly5NCV{0,   0,   0,    0,    0,    0,    0,    0,    0,     0,    0,
                                    0,   0,   0,    0,    0,    1,    1,    1,    1,     1,    1,
                                    1,   1,   -1e8, -1e8, -1e8, -1e8, -1e8, -1e8, -1000, -1e8, 1e8,
                                    1e8, 1e8, 1e8,  1e8,  1e8,  1000, 1e8};  // poly 5
// priors for adv dir down
std::vector<double> prRFreqPoly1Down{0, 0,    0,    0,    0,   0, 1, 1,
                                     1, -1e8, -1e8, -1e8, 1e8, 0, 0};  // poly 1
std::vector<double> prRFreqPoly2Down{0, 0, 0,    0,    0,    0,    0,   0, 1, 1,
                                     1, 1, -1e8, -1e8, -1e8, -1e8, 1e8, 0, 0, 0};  // poly 2
std::vector<double> prRFreqPoly3Down{
    0, 0, 0,    0,    0,    0,    0,    0,   0, 0, 1, 1, 1,
    1, 1, -1e8, -1e8, -1e8, -1e8, -1e8, 1e8, 0, 0, 0, 0
};  // poly 3
std::vector<double> prRFreqPoly4Down{0,    0,    0,    0,    0,   0, 0, 0, 0,    0,
                                     0,    0,    1,    1,    1,   1, 1, 1, -1e8, -1e8,
                                     -1e8, -1e8, -1e8, -1e8, 1e8, 0, 0, 0, 0,    0};  // poly 4
std::vector<double> prRFreqPoly5Down{0,    0,    0,    0,    0,   0, 0, 0, 0, 0,    0,    0,
                                     0,    0,    1,    1,    1,   1, 1, 1, 1, -1e8, -1e8, -1e8,
                                     -1e8, -1e8, -1e8, -1e8, 1e8, 0, 0, 0, 0, 0,    0};  // poly 5
std::vector<double> prRFreqPoly1NCVDown{0, 0, 0,    0,    0,    0,    0,   0, 1,    1,
                                        1, 1, -1e8, -1e8, 1000, -1e8, 1e8, 0, 1000, 0};  // poly 1
std::vector<double> prRFreqPoly2NCVDown{
    0, 0, 0,    0,    0,    0,     0,    0,   0, 0, 1,    1, 1,
    1, 1, -1e8, -1e8, -1e8, -1000, -1e8, 1e8, 0, 0, 1000, 0
};  // poly 2
std::vector<double> prRFreqPoly3NCVDown{0,    0,    0,     0,    0,   0, 0, 0, 0,    0,
                                        0,    0,    1,     1,    1,   1, 1, 1, -1e8, -1e8,
                                        -1e8, -1e8, -1000, -1e8, 1e8, 0, 0, 0, 1000, 0};  // poly 3
std::vector<double> prRFreqPoly4NCVDown{
    0, 0, 0, 0,    0,    0,    0,    0,    0,     0,    0,   0, 0, 0, 1, 1,    1, 1,
    1, 1, 1, -1e8, -1e8, -1e8, -1e8, -1e8, -1000, -1e8, 1e8, 0, 0, 0, 0, 1000, 0
};  // poly 4
std::vector<double> prRFreqPoly5NCVDown{0, 0, 0,    0,    0,    0,    0,    0,    0,     0,    0,
                                        0, 0, 0,    0,    0,    1,    1,    1,    1,     1,    1,
                                        1, 1, -1e8, -1e8, -1e8, -1e8, -1e8, -1e8, -1000, -1e8, 1e8,
                                        0, 0, 0,    0,    0,    1000, 0};  // poly 5

// priors for adv dir up
std::vector<double> prRFreqPoly1Up{0, 0, 0, 0, 0, 0, 1, 1, 1, -1e8, 0, 0, 1e8, 1e8, 1e8};  // poly 1
std::vector<double> prRFreqPoly2Up{0, 0, 0,    0, 0, 0, 0,   0,   1,   1,
                                   1, 1, -1e8, 0, 0, 0, 1e8, 1e8, 1e8, 1e8};  // poly 2
std::vector<double> prRFreqPoly3Up{0, 0, 0,    0, 0, 0, 0, 0,   0,   0,   1,   1,  1,
                                   1, 1, -1e8, 0, 0, 0, 0, 1e8, 1e8, 1e8, 1e8, 1e8};  // poly 3
std::vector<double> prRFreqPoly4Up{
    0, 0, 0, 0,    0, 0, 0, 0, 0, 0,   0,   0,   1,   1,   1,
    1, 1, 1, -1e8, 0, 0, 0, 0, 0, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8
};  // poly 4
std::vector<double> prRFreqPoly5Up{0, 0, 0, 0, 0,   0,   0,   0,   0,   0,    0,  0,
                                   0, 0, 1, 1, 1,   1,   1,   1,   1,   -1e8, 0,  0,
                                   0, 0, 0, 0, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8,  1e8};  // poly 5

// old
std::vector<double> prRFreqPoly1NCVUp{0, 0, 0,    0, 0,     0, 0,   0,   1,    1,
                                      1, 1, -1e8, 0, -1000, 0, 1e8, 1e8, 1000, 1e8};  // poly 1
// new
std::vector<double> prRFreqPoly2NCVUp{
    0, 0, 0,    0, 0, 0,     0,     0,   0,   0,   1,    1,   1,
    1, 1, -1e8, 0, 0, -1000, -1000, 1e8, 1e8, 1e8, 1000, 1000
};  // poly 2
std::vector<double> prRFreqPoly3NCVUp{0,     0,     0,   0,   0,   0,   0,    0,    0, 0, 0,
                                      0,     1,     1,   1,   1,   1,   1,    -1e8, 0, 0, 0,
                                      -1000, -1000, 1e8, 1e8, 1e8, 1e8, 1000, 1000};  // poly 3
std::vector<double> prRFreqPoly4NCVUp{
    0, 0, 0, 0,    0, 0, 0, 0, 0,     0,     0,   0,   0,   0,   1,   1,    1,   1,
    1, 1, 1, -1e8, 0, 0, 0, 0, -1000, -1000, 1e8, 1e8, 1e8, 1e8, 1e8, 1000, 1000
};  // poly 4
std::vector<double> prRFreqPoly5NCVUp{0,   0,   0,    0,   0,   0,    0,   0, 0,     0,     0,
                                      0,   0,   0,    0,   0,   1,    1,   1, 1,     1,     1,
                                      1,   1,   -1e8, 0,   0,   0,    0,   0, -1000, -1000, 1e8,
                                      1e8, 1e8, 1e8,  1e8, 1e8, 1000, 1000};  // poly 5

// MS cancer
double prRFreqMultistageCancerG[] = {0, -17, 0, -18, 18};
double prRFreqMultistageCancerB[] = {0, 0.1, 0, 0, 1e4};
