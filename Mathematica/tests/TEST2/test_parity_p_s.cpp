#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include "../.././AxolotE_Mat.h"

bool run_test_l1_l0() {
    const double alpha = 1.62;
    const double beta  = 0.88;

    // Coordinates A and B
    const std::array<double,3> A_orig = {  0.45, -1.20,  0.85 };
    const std::array<double,3> B_orig = { -1.35,  0.70, -0.40 };

    // Inverted coordinates A' = -A, B' = -B
    const std::array<double,3> A_inv  = { -A_orig[0], -A_orig[1], -A_orig[2] };
    const std::array<double,3> B_inv  = { -B_orig[0], -B_orig[1], -B_orig[2] };

    // Expansion coefficients
    std::vector<double> E_orig = ECoeff1And0(alpha, beta, A_orig, B_orig);
    std::vector<double> E_inv  = ECoeff1And0(alpha, beta, A_inv,  B_inv);

    if (E_orig.empty() || E_orig.size() != E_inv.size()) {
        std::cerr << "[FAIL] Vector size mismatch or empty result for l=" 
                  << 1 << ", lp=" << 0 << "\n";
        return false;
    }

    // Verify magnitude invariance: |E_orig[i]| == |E_inv[i]|
    for (size_t i = 0; i < E_orig.size(); ++i) {
        double mag_orig = std::abs(E_orig[i]);
        double mag_inv  = std::abs(E_inv[i]);

        if (std::abs(mag_orig - mag_inv) > 1e-12) {
            std::cerr << "[FAIL] Parity magnitude broken for l=" << 1 << ", lp=" << 0 
                      << " at index " << i 
                      << " (|orig|=" << mag_orig << ", |inv|=" << mag_inv << ")\n";
            return false;
        }
    }

    return true;
}
