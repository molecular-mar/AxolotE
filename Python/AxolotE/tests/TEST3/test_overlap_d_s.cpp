#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include "../.././AxolotE_Py.h"

// Forward declaration for targeted generated function
std::vector<double> ECoeff2And0(double expAlpha, double expBeta,
                                      std::array<double,3> muCoords, 
                                      std::array<double,3> nuCoords);

bool run_test_l2_l0() {
    const double alpha = 1.25;
    const double beta  = 0.75;

    // Single center position (A = B)
    const std::array<double,3> center = { 0.5, -1.2, 2.3 };

    // Evaluate generated expansion coefficients at coincident centers
    std::vector<double> E_coeffs = ECoeff2And0(alpha, beta, center, center);

    const int l1 = 2;
    const int l2 = 0;

    const int num_m1 = 2 * l1 + 1;
    const int num_m2 = 2 * l2 + 1;
    const int num_pairs = num_m1 * num_m2;

    if (E_coeffs.empty() || E_coeffs.size() % num_pairs != 0) {
        std::cerr << "[FAIL] Invalid coefficient vector size for l=" 
                  << l1 << ", lp=" << l2 << "\n";
        return false;
    }

    const size_t n_hermite = E_coeffs.size() / num_pairs;

    // Verify single-center orthogonality and diagonal overlap
    for (int m1_idx = 0; m1_idx < num_m1; ++m1_idx) {
        for (int m2_idx = 0; m2_idx < num_m2; ++m2_idx) {
            int pair_idx = m1_idx * num_m2 + m2_idx;
            size_t e000_idx = pair_idx * n_hermite;

            double e000_val = E_coeffs[e000_idx];

            if (l1 != l2) {
                // Different angular momentum shells (l1 != l2): E_{0,0,0} must be 0
                if (std::abs(e000_val) > 1e-12) {
                    std::cerr << "[FAIL] Cross-shell single-center orthogonality violated for l1=" 
                              << l1 << ", l2=" << l2 << " at pair (" 
                              << m1_idx << "," << m2_idx << "): E000 = " 
                              << e000_val << "\n";
                    return false;
                }
            } else {
                // Same angular momentum shell (l1 == l2)
                if (m1_idx != m2_idx) {
                    // Orthogonal m components: E_{0,0,0} must be 0
                    if (std::abs(e000_val) > 1e-12) {
                        std::cerr << "[FAIL] Single-center m-component orthogonality violated for l=" 
                                  << l1 << " at pair (" << m1_idx << "," << m2_idx 
                                  << "): E000 = " << e000_val << "\n";
                        return false;
                    }
                } else {
                    // Identical orbital (m1 == m2): E_{0,0,0} must be strictly positive
                    if (e000_val <= 1e-12) {
                        std::cerr << "[FAIL] Single-center diagonal overlap non-positive for l=" 
                                  << l1 << " at component " << m1_idx 
                                  << ": E000 = " << e000_val << "\n";
                        return false;
                    }
                }
            }
        }
    }

    std::cout << "[PASS] Shell pair (" << 2 << ", " << 0 
              << ") evaluated successfully.\n";

    return true;
}
