// test_main.cpp

#include "test_main.h"

#include "integration_tests.h"
#include "unit_tests.h"

int main() {
  run_all_unitTests();
  run_all_integrationTests();
}
