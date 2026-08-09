# AxolotE: Python implementation

Files for the generation of code with the E Coefficients expressions.

---

## Directory Contents

* **`ERecurrence.py`**: Core algorithm containing the fundamental McMurchie-Davidson recurrence relations for $E$ coefficients.
* **`EGenerator.py`**: Symbolic generator module containing instructions for formatting and writing explicit coefficient expressions (e.g., in Mathematica-compatible format).
* **`CppCodeGenerator.py`**: The primary Python script that executes the automated code generation, converting mathematical recurrences into explicit C++ functions.
* **`starters/`**: Directory containing C++ template files used by `CppCodeGenerator.py` during code synthesis.
* **`tests/`**: Unit testing suite directory containing validation suites for global translational invariance (`TEST1`), coordinate inversion parity symmetry (`TEST2`), and single-center overlap/orthogonality (`TEST3`).

---

# Usage Instructions

### 1. Requirements & Prerequisites
* Python 3.9 or higher
* C++11 (or higher) compliant compiler 
* SymPy 1.2 or higher

---

### 2. Generating C++ Source Code

To generate or update the explicit C++ routines, first modify `CppCodeGenerator.py` to indicate the ($\ell,\ell'$) values to be considered, following the comments inside the file. Then, execute `CppCodeGenerator.py` from the command line:

```bash
python CppCodeGenerator.py
```

After correct execution, files `AxolotE_Py.h` and `AxolotE_Py.cpp` will be generated ready to be used.
---

### 3. Integrating Generated Code into Your Project

Include `AxolotE_Py.h` in your project source files and compile `AxolotE_Py.cpp` alongside your application:

```cpp
#include "AxolotE_Py.h"

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
g++ -std=c++11 main.cpp AxolotE_Py.cpp -o my_integral_program
```

---

### 4. Running Tests

To verify the spatial and physical consistency of the generated C++ routines, as well as to find usage examples, navigate to the desired test directory under `tests/` and follow the instructions in its respective `README.md`
