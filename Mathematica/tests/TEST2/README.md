# Coordinate Inversion Parity Symmetry Unit Tests

This directory provides automated unit tests to verify the spatial inversion parity symmetry of code-generated expansion coefficients $E$.

---

### Theory 

The $E$ coefficients depend strictly on relative displacement vectors:

$$\mathbf{D_{PA}} = \mathbf{R_P} - \mathbf{R_A} = \frac{\beta}{\alpha + \beta}(\mathbf{B} - \mathbf{A})$$

$$\mathbf{D_{PB}} = \mathbf{R_P} - \mathbf{R_B} = -\frac{\alpha}{\alpha + \beta}(\mathbf{B} - \mathbf{A})$$

Under a coordinate inversion transformation about the origin ($\mathbf{A}' = -\mathbf{A}$ and $\mathbf{B}' = -\mathbf{B}$), the relative displacement vectors invert their direction ($\mathbf{D_{PA}}' = -\mathbf{D_{PA}}$ and $\mathbf{D_{PB}}' = -\mathbf{D_{PB}}$). 

Because every polynomial term contributing to a given $E_{t,u,v}$ coefficient possesses a well-defined parity under coordinate inversion, inversion alters at most the overall sign of any coefficient. Consequently, the absolute magnitude of every calculated Hermite expansion coefficient must remain strictly invariant:

$$\left| E_{t,u,v}^{\ell_1 m_1, \ell_2 m_2}(-\mathbf{A}, -\mathbf{B}, \alpha, \beta) \right| = \left| E_{t,u,v}^{\ell_1 m_1, \ell_2 m_2}(\mathbf{A}, \mathbf{B}, \alpha, \beta) \right|$$

Testing this magnitude symmetry element-by-element across all generated $E$-coefficient vectors guarantees that the automated code generator correctly handles term sign assignments and coordinate subtraction ordering.

---

### Folder Contents

Two files are provided to generate the source code for the tests:

* `parity_template.cpp.in`: Reference C++ unit test template. It defines test geometry and inverted vectors, and contains an `@AXOLOTE_DIR@` placeholder for the generated code directory, and `@L1@` and `@L2@` placeholders for shell angular momentum quantum numbers.
* `generate_parity_tests.sh`: Bash generation script. It substitutes `@AXOLOTE_DIR@` for the generated code directory, and `@L1@` and `@L2@` across all combination pairs of $s, p, d$ shells ($\ell, \ell' \in \{0, 1, 2\}$) to produce individual C++ unit test files and a main test suite runner `suite_parity_test.cpp`.

Pre-generated tests for $s, p,$ and $d$ shells are included:

* `test_parity_*.cpp`: C++ functions to perform the test for each individual $\ell,\ell'$ combination.
* `suite_parity_test.cpp`: C++ main function to execute all the tests.

---

### Instructions

#### 1. Compilation

Compile the generated test suite files together with the generated coefficient source code:

```bash
g++ -std=c++11 \
    -I../.. \
    suite_parity_test.cpp \
    test_parity_*.cpp \
    ../../AxolotE_Mat.cpp \
    -o run_parity_tests
```

#### 2. Execution

Execute the test binary to evaluate all available tests:

```bash
./run_parity_tests
```

#### 3. Generate Test Suite

Modify as desired the `generate_parity_tests.sh` script to generate tests for different $(\ell,\ell')$ combinations.  

Run the shell script to create the test files and master runner inside the directory:

```bash
bash generate_parity_tests.sh
```