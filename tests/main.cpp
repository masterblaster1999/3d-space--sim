#include <iostream>

int test_orbit();
int test_streaming();
int test_proc();
int test_economy();
int test_route_planner();
int test_savegame();
int test_missions();
int test_nav();

int main() {
  int fails = 0;

  fails += test_orbit();
  fails += test_streaming();
  fails += test_proc();
  fails += test_economy();
  fails += test_route_planner();
  fails += test_savegame();
  fails += test_missions();
  fails += test_nav();

  if (fails == 0) {
    std::cout << "[stellar_tests] ALL PASS\n";
    return 0;
  }

  std::cerr << "[stellar_tests] FAILS=" << fails << "\n";
  return 1;
}
