# Single-Center Overlap and Orthogonality Unit Tests

This directory provides automated unit tests to verify single-center orbital orthogonality and overlap properties of code-generated expansion coefficients $E$.

---

### Theory 

In McMurchie-Davidson theory, the spatial overlap integral $S_{\mu\nu}$ between two basis functions $\chi_\mu$ and $\chi_\nu$ reduces to the zero-order Hermite coefficient ($t=0, u=0, v=0$):

$$S_{\mu\nu} = \int \chi_\mu(\mathbf{r}) \chi_\nu(\mathbf{r}) \, d\mathbf{r} = E_{0,0,0}^{\ell_1 m_1, \ell_2 m_2}(\mathbf{A}, \mathbf{B}, \alpha, \beta) \cdot \left(\frac{\pi}{\alpha + \beta}\right)^{3/2}$$

When both Gaussian functions share an identical center ($\mathbf{A} = \mathbf{B}$), real solid harmonic functions $r^{\ell} S_{\ell m}(\hat{r})$ obey exact analytical orthogonality relations over angular integration:

$$E_{0,0,0}^{\ell_1 m_1, \ell_2 m_2}(\mathbf{A}, \mathbf{A}, \alpha, \beta) \propto \delta_{\ell_1 \ell_2} \, \delta_{m_1 m_2}$$

Consequently, evaluating coefficients on a single center imposes two fundamental conditions:

1. **Orthogonality:** $E_{0,0,0}^{\ell_1 m_1, \ell_2 m_2} = 0$ whenever $\ell_1 \neq \ell_2$ or $m_1 \neq m_2$.
2. **Positivity:** $E_{0,0,0}^{\ell_1 m_1, \ell_1 m_1} > 0$ for identical basis functions ($\ell_1 = \ell_2$ and $m_1 = m_2$).

Verifying these selection rules guarantees that the automated code generator correctly preserves the underlying angular momentum coupling and real solid harmonic orthogonality.

---

### Folder Contents

Two files are provided to generate the source code for the tests:

* `overlap_template.cpp.in`: Reference C++ unit test template. It evaluates coefficients at coincident nuclear centers ($\mathbf{A} = \mathbf{B}$), and contains an `@AXOLOTE_DIR@` placeholder for the generated code directory, and `@L1@` and `@L2@` placeholders for shell angular momentum quantum numbers.
* `generate_overlap_tests.sh`: Bash generation script. It substitutes `@AXOLOTE_DIR@` for the generated code directory, and `@L1@` and `@L2@` across all combination pairs of $s, p, d$ shells ($\ell, \ell' \in \{0, 1, 2\}$) to produce individual C++ unit test files and a main test suite runner `suite_overlap_test.cpp`.

Pre-generated tests for $s, p,$ and $d$ shells are included:

* `test_overlap_*.cpp`: C++ functions to perform the test for each individual $\ell,\ell'$ combination.
* `suite_overlap_test.cpp`: C++ main function to execute all the tests.

---

### Instructions

#### 1. Compilation

Compile the generated test suite files together with the generated coefficient source code:

```bash
g++ -std=c++11 \
    -I../.. \
    suite_overlap_test.cpp \
    test_overlap_*.cpp \
    ../../AxolotE_Mat.cpp \
    -o run_overlap_tests 
```

#### 2. Execution

Execute the test binary to evaluate all available tests:

```bash
./run_overlap_tests
```

#### 3. Generate Test Suite

Modify as desired the `generate_overlap_tests.sh` script to generate tests for different $(\ell,\ell')$ combinations.  

Run the shell script to create the test files and master runner inside the directory:

```bash
bash generate_overlap_tests.sh
```
