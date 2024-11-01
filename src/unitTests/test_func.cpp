#include "test_func.h"

bool essentiallyEqual(double a, double b, double epsilon) {
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}
