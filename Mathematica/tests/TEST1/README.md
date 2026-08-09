# Global Translational Invariance Unit Tests

This directory provides automated unit tests to verify the spatial symmetry of code-generated expansion coefficients $E$.

---

### Theory 

The E coefficients depend strictly on relative displacement vectors:


$$\mathbf{D_PA} = \mathbf{R_P} - \mathbf{R_A} = \frac{\beta}{\alpha + \beta}(\mathbf{B} - \mathbf{A})$$

$$\mathbf{D_PB} = \mathbf{R_P} - \mathbf{R_B} = -\frac{\alpha}{\alpha + \beta}(\mathbf{B} - \mathbf{A})$$

Under a rigid global translation by an arbitrary vector $\mathbf{r}_0$ ($\mathbf{A}' = \mathbf{A} + \mathbf{r}_0$ and $\mathbf{B}' = \mathbf{B} + \mathbf{r}_0$), the inter-center distance remains unchanged:


$$\mathbf{B}' - \mathbf{A}' = (\mathbf{B} + \mathbf{r}_0) - (\mathbf{A} + \mathbf{r}_0) = \mathbf{B} - \mathbf{A}$$

Consequently, every calculated Hermite expansion coefficient must remain strictly invariant:

$$E_{t,u,v}^{\ell_1 m_1, \ell_2 m_2}(\mathbf{A} + \mathbf{r}_0, \mathbf{B} + \mathbf{r}_0, \alpha, \beta) = E_{t,u,v}^{\ell_1 m_1, \ell_2 m_2}(\mathbf{A}, \mathbf{B}, \alpha, \beta)$$

Testing this symmetry element-by-element across all generated $E$-coefficient vectors guarantees that the automated code generator correctly unrolled distance expressions without hardcoding absolute coordinates or introducing sign errors.

---

### Folder Contents

Two files are provided to generate the source code for the tests:

* `translation_template.cpp.in`: Reference C++ unit test template. It defines test geometry and displacement vectors, and contains a `@AXOLOTE_DIR@` placeholder for the generated code directory, and `@L1@` and `@L2@` placeholders for shell angular momentum quantum numbers.
* `generate_translation_tests.sh`: Bash generation script. It substitutes `@AXOLOTE_DIR@` for the generated code directory, and `@L1@` and `@L2@` across all combination pairs of $s, p, d$ shells ($\ell, \ell' \in \{0, 1, 2\}$) to produce individual C++ unit test files and a main test suite runner `suite_translation_test.cpp`.

Pre-generated tests for $s, p,$ and $d$ shells are included:

* `test_translation_*.cpp`: C++ functions to perform the test for each individual $\ell,\ell'$ combination.
* `suite_translation_test.cpp`: C++ main function to execute all the tests.

---

### Instructions

#### 1. Compilation

Compile the generated test suite files together with the generated coefficient source code:

```bash
g++ -std=c++11 \
    -I../.. \
    suite_translation_test.cpp \
    test_translation_*.cpp \
    ../../AxolotE_Mat.cpp \
    -o run_translation_tests 
```

#### 2. Execution

Execute the test binary to evaluate all available tests:

```bash
./run_translation_tests

```

#### 3. Generate Test Suite

Modify as desired the `generate_translation_tests.sh` script to generate tests for different $(\ell,\ell')$ combinations.  

Run the shell script to create the test files and master runner inside the directory:

```bash
bash generate_translaiton_tests.sh
```


