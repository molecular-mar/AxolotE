#!/usr/bin/env bash
#Script to generate unit tests for the translation invariance test.
#Adapt in order to generate tests for different l and lp combinations

#Exit if any error is generated at any point
set -e

AXOLOTE_DIR="../../."
TEMPLATE_FILE="translation_template.cpp.in"
OUTPUT_DIR="."
MAIN_RUNNER="${OUTPUT_DIR}/suite_translation_test.cpp"

if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "Error: Template file '$TEMPLATE_FILE' not found."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Shell type mapping for better output
declare -A SHELL_NAMES=([0]="s" [1]="p" [2]="d")

# Generation of full tests program
cat << 'EOF' > "$MAIN_RUNNER"
#include <iostream>
#include <cstdlib>

EOF

for L1 in 0 1 2; do
    for L2 in 0 1 2; do
        echo "bool run_test_l${L1}_l${L2}();" >> "$MAIN_RUNNER"
    done
done

cat << 'EOF' >> "$MAIN_RUNNER"

int main() {
    std::cout << "========================================\n";
    std::cout << " Running Global Invariance Unit Tests   \n";
    std::cout << "========================================\n";

    int passed = 0;
    int total = 0;

EOF

# Generate individual test C++ files for all s, p, d combinations (0, 1, 2)
for L1 in 0 1 2; do
    for L2 in 0 1 2; do
        S1=${SHELL_NAMES[$L1]}
        S2=${SHELL_NAMES[$L2]}
        TEST_FILE="${OUTPUT_DIR}/test_translation_${S1}_${S2}.cpp"

        # Replace placeholders in template to build individual test file
        sed -e "s#@AXOLOTE_DIR@#${AXOLOTE_DIR}#g" -e "s/@L1@/${L1}/g" -e "s/@L2@/${L2}/g" "$TEMPLATE_FILE" > "$TEST_FILE"
        echo "Generated: ${TEST_FILE} (${S1}-${S2} pair)"

        # Append execution call to main runner
        cat << EOF >> "$MAIN_RUNNER"
    total++;
    if (run_test_l${L1}_l${L2}()) passed++;
EOF
    done
done

cat << 'EOF' >> "$MAIN_RUNNER"

    std::cout << "========================================\n";
    std::cout << " Summary: " << passed << "/" << total << " tests passed.\n";
    std::cout << "========================================\n";

    return (passed == total) ? EXIT_SUCCESS : EXIT_FAILURE;
}
EOF

echo ""
echo "Generation complete. All test files generated in '${OUTPUT_DIR}/'."
