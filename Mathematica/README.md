# AxolotE: Mathematica implementation

File *AxolotE.m* is a Mathematica Script with all the tools requiered for generating and testing the E Coefficients. Detailed information can be found inside the Notebook.

---

## Directory Contents

* **`AxolotE.m`**: Mathematica Script with the instructions and commands to generate the C++ code and to validate the generated algebraic expressions.
* **`AxolotE_Mat.cpp` and  `AxolotE_Mat.h`**: Pre-generated C++ code for evaluation of $E$ coefficients, up to $(\ell,\ell')=(2,2)$.
* **`tests/`**: Unit testing suite directory containing validation suites for global translational invariance (`TEST1`), coordinate inversion parity symmetry (`TEST2`), and single-center overlap/orthogonality (`TEST3`).

---

# Usage Instructions

### 1. Requirements & Prerequisites
* Wolfram Mathematica (tested with version 12.0)
* C++11 (or higher) compliant compiler 

### 2. Generating C++ Source Code

To generate or update the explicit C++ routines, modify and run the `AxolotE.m` file.

After correct execution, files `AxolotE_Mat.h` and `AxolotE_Mat.cpp` will be generated ready to be used.
---

### 3. Integrating Generated Code into Your Project

Include `AxolotE_Mat.h` in your project source files and compile `AxolotE_Mat.cpp` alongside your application:

```cpp
#include "AxolotE_Mat.h"

int main() {
    // Example: Evaluate (s, p) pair coefficients
    double alpha = 1.2, beta = 0.8;
    std::array<double, 3> A = {0.0, 0.0, 0.0};
    std::array<double, 3> B = {1.0, 0.5, -0.2};

    std::vector<double> E = ECoeff0And1(alpha, beta, A, B);

    return 0;
}
```

Compile with your target program:

```bash
g++ -std=c++11 main.cpp AxolotE_Mat.cpp -o my_integral_program
```

---

### 4. Running Tests

To verify the spatial and physical consistency of the generated C++ routines, as well as to find usage examples, navigate to the desired test directory under `tests/` and follow the instructions in its respective `README.md`
