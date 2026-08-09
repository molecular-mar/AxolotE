#include <iostream>
#include <vector>
#include <chrono>
#include <array>
#include <iomanip>
#include <cmath>

// 1. Pre-include standard library headers so their header guards execute in global scope.

// 2. Wrap the Python implementation in the PyImpl namespace
namespace PyImpl {
    #include "../Python/AxolotE/AxolotE_Py.h"
    #include "../Python/AxolotE/AxolotE_Py.cpp"
}

// 3. Wrap the Mathematica implementation in the MatImpl namespace
namespace MatImpl {
    #include "../Mathematica/AxolotE_Mat.h"
    #include "../Mathematica/AxolotE_Mat.cpp"
}

using Coords = std::array<double, 3>;
using ECoeffFunc = std::vector<double>(*)(double, double, Coords, Coords);

// =========================================================================
// Benchmark Harness
// =========================================================================
double benchmarkImplementation(ECoeffFunc func, 
                               double alpha, double beta, 
                               const Coords& ACoords, const Coords& BCoords, 
                               int iterations = 10000) {
    if (!func) return 0.0;

    // Warm-up
    for (int i = 0; i < 100; ++i) {
        auto dummy = func(alpha, beta, ACoords, BCoords);
        volatile double block_opt = dummy.empty() ? 0.0 : dummy[0]; 
    }

    // Timing loop
    volatile double sink = 0.0; 
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        auto result = func(alpha, beta, ACoords, BCoords);
        if (!result.empty()) {
            sink = sink + result[0]; 
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::micro> duration = end - start;
    return duration.count() / iterations;
}

// =========================================================================
// Main Benchmark Routine
// =========================================================================
int main() {
    double alpha{0.53318702000};
    double beta{0.7865398200};
    Coords ACoords{0.724, 1.255, 7.7};
    Coords BCoords{1.449, 2.509, 7.3};
    int iterations = 1000;

    // Approach 1: Python Implementation
    ECoeffFunc approach1[5][5] = {nullptr};
    approach1[0][0] = &PyImpl::ECoeff0And0;
    approach1[0][1] = &PyImpl::ECoeff0And1;
    approach1[0][2] = &PyImpl::ECoeff0And2;
    //approach1[0][3] = &PyImpl::ECoeff0And3;
    //approach1[0][4] = &PyImpl::ECoeff0And4;
    approach1[1][0] = &PyImpl::ECoeff1And0;
    approach1[1][1] = &PyImpl::ECoeff1And1;
    approach1[1][2] = &PyImpl::ECoeff1And2;
    //approach1[1][3] = &PyImpl::ECoeff1And3;
    //approach1[1][4] = &PyImpl::ECoeff1And4;
    approach1[2][0] = &PyImpl::ECoeff2And0;
    approach1[2][1] = &PyImpl::ECoeff2And1;
    approach1[2][2] = &PyImpl::ECoeff2And2;
    //approach1[2][3] = &PyImpl::ECoeff2And3;
    //approach1[2][4] = &PyImpl::ECoeff2And4;
    /*approach1[3][0] = &PyImpl::ECoeff3And0;
    approach1[3][1] = &PyImpl::ECoeff3And1;
    approach1[3][2] = &PyImpl::ECoeff3And2;
    approach1[3][3] = &PyImpl::ECoeff3And3;
    approach1[3][4] = &PyImpl::ECoeff3And4;
    approach1[4][0] = &PyImpl::ECoeff4And0;
    approach1[4][1] = &PyImpl::ECoeff4And1;
    approach1[4][2] = &PyImpl::ECoeff4And2;
    approach1[4][3] = &PyImpl::ECoeff4And3;
    approach1[4][4] = &PyImpl::ECoeff4And4;*/

    // Approach 2: Mathematica Implementation
    ECoeffFunc approach2[5][5] = {nullptr};
    approach2[0][0] = &MatImpl::ECoeff0And0;
    approach2[0][1] = &MatImpl::ECoeff0And1;
    approach2[0][2] = &MatImpl::ECoeff0And2;
    //approach2[0][3] = &MatImpl::ECoeff0And3;
    //approach2[0][4] = &MatImpl::ECoeff0And4;
    approach2[1][0] = &MatImpl::ECoeff1And0;
    approach2[1][1] = &MatImpl::ECoeff1And1;
    approach2[1][2] = &MatImpl::ECoeff1And2;
    //approach2[1][3] = &MatImpl::ECoeff1And3;
    //approach2[1][4] = &MatImpl::ECoeff1And4;
    approach2[2][0] = &MatImpl::ECoeff2And0;
    approach2[2][1] = &MatImpl::ECoeff2And1;
    approach2[2][2] = &MatImpl::ECoeff2And2;
    //approach2[2][3] = &MatImpl::ECoeff2And3;
    //approach2[2][4] = &MatImpl::ECoeff2And4;
    /*approach2[3][0] = &MatImpl::ECoeff3And0;
    approach2[3][1] = &MatImpl::ECoeff3And1;
    approach2[3][2] = &MatImpl::ECoeff3And2;
    approach2[3][3] = &MatImpl::ECoeff3And3;
    approach2[3][4] = &MatImpl::ECoeff3And4;
    approach2[4][0] = &MatImpl::ECoeff4And0;
    approach2[4][1] = &MatImpl::ECoeff4And1;
    approach2[4][2] = &MatImpl::ECoeff4And2;
    approach2[4][3] = &MatImpl::ECoeff4And3;
    approach2[4][4] = &MatImpl::ECoeff4And4;*/

    std::cout << "Benchmarking AxolotE implementations (Times in microseconds)\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << " L1 | L2 | Python (us)     | Mathematica (us)| Diff (%)\n";
    std::cout << "--------------------------------------------------------\n";

    for (int l1 = 0; l1 <= 2; ++l1) {
        for (int l2 = 0; l2 <= 2; ++l2) {
            if (!approach1[l1][l2] || !approach2[l1][l2]) continue;

            double time1 = benchmarkImplementation(approach1[l1][l2], alpha, beta, ACoords, BCoords, iterations);
            double time2 = benchmarkImplementation(approach2[l1][l2], alpha, beta, ACoords, BCoords, iterations);
            
            double diff = ((time2 - time1) / time2) * 100.0;

            std::cout << "  " << l1 << " |  " << l2 << " | "
                      << std::fixed << std::setprecision(4) << std::setw(15) << time1 << " | "
                      << std::setw(15) << time2 << " | "
                      << std::setw(7) << diff << "%\n";
        }
    }

    return 0;
}
