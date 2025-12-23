#include <iostream>

int test_orbit();
int test_proc();

int main() {
  int failures = 0;

  failures += test_orbit();
  failures += test_proc();

  if (failures == 0) {
    std::cout << "[stellar_tests] All tests passed.\n";
    return 0;
  }

  std::cerr << "[stellar_tests] Failures: " << failures << "\n";
  return 1;
}
