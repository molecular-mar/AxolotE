#include <iostream>
#include <cstdlib>

bool run_test_l0_l0();
bool run_test_l0_l1();
bool run_test_l0_l2();
bool run_test_l1_l0();
bool run_test_l1_l1();
bool run_test_l1_l2();
bool run_test_l2_l0();
bool run_test_l2_l1();
bool run_test_l2_l2();

int main() {
    std::cout << "========================================\n";
    std::cout << " Running Single-Center Overlap Tests   \n";
    std::cout << "========================================\n";

    int passed = 0;
    int total = 0;

    total++;
    if (run_test_l0_l0()) passed++;
    total++;
    if (run_test_l0_l1()) passed++;
    total++;
    if (run_test_l0_l2()) passed++;
    total++;
    if (run_test_l1_l0()) passed++;
    total++;
    if (run_test_l1_l1()) passed++;
    total++;
    if (run_test_l1_l2()) passed++;
    total++;
    if (run_test_l2_l0()) passed++;
    total++;
    if (run_test_l2_l1()) passed++;
    total++;
    if (run_test_l2_l2()) passed++;

    std::cout << "========================================\n";
    std::cout << " Summary: " << passed << "/" << total << " tests passed.\n";
    std::cout << "========================================\n";

    return (passed == total) ? EXIT_SUCCESS : EXIT_FAILURE;
}
