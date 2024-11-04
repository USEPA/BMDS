#include <cmath>

#define expect_true(arg)                                                                      \
  do {                                                                                        \
    if (!(arg)) {                                                                             \
      std::cout << "Unexpected false at " << __FILE__ << ", " << __LINE__ << ", " << __func__ \
                << ": " << #arg << std::endl;                                                 \
    }                                                                                         \
  } while (false);

bool essentiallyEqual(double a, double b, double epsilon);
