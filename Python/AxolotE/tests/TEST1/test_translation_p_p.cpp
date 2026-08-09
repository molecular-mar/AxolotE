#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include "../.././AxolotE_Py.h"

bool run_test_l1_l1() {
    const double alpha = 1.45;
    const double beta  = 0.75;

    const std::array<double,3> mu_orig  = { 0.2, -1.1,  0.5 };
    const std::array<double,3> nu_orig  = { 1.8,  0.3, -0.4 };
    const std::array<double,3> shift    = { 12.5, -45.0, 102.75 };

    const std::array<double,3> mu_shift = { mu_orig[0] + shift[0], mu_orig[1] + shift[1], mu_orig[2] + shift[2] };
    const std::array<double,3> nu_shift = { nu_orig[0] + shift[0], nu_orig[1] + shift[1], nu_orig[2] + shift[2] };

    std::vector<double> E_orig  = ECoeff1And1(alpha, beta, mu_orig, nu_orig);
    std::vector<double> E_shift = ECoeff1And1(alpha, beta, mu_shift, nu_shift);

    if (E_orig.empty() || E_orig.size() != E_shift.size()) {
        std::cerr << "[FAIL] Size mismatch or empty vector for l=" << 1 << ", lp=" << 1 << "\n";
        return false;
    }

    for (size_t i = 0; i < E_orig.size(); ++i) {
        if (std::abs(E_orig[i] - E_shift[i]) > 1e-12) {
            std::cerr << "[FAIL] Invariance broken for l=" << 1 << ", lp=" << 1 
                      << " at index " << i << "\n";
            return false;
        }
    }

    std::cout << "[PASS] Shell Pair (" << 1 << ", " << 1 << ") | Tested " 
              << E_orig.size() << " coefficients.\n";
    return true;
}
