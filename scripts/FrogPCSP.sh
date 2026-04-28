#!/usr/bin/env bash
# ==============================
# Pipeline runner for cleavage site prediction
# Runs three Python scripts in sequence:
#   1. svm-predictor.py <input.fasta>
#   2. make-fragments.py <input.fasta> results/predictions.csv
#   3. compute-physicochemical-properties.py results/fragments.fasta
#
# Usage:
#   ./FrogPCSP.sh <input.fasta>
# ==============================

# Check input argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.fasta>"
    exit 1
fi

INPUT_FASTA="$1"

# Paths to Python scripts
PREDICTOR="scripts/svm-predictor.py"
FRAGMENTS="scripts/make-fragments.py"
PHYCH="scripts/compute-physicochemical-properties.py"
SCAT="scripts/scatterplot.py"

# Output files
OUTPUT_CSV="results/predictions.csv"
OUTPUT_FRAG="results/fragments.fasta"
OUTPUT_PHYSCH="results/phychem.tsv"
OUTPUT_SCAT="results/scatterplot.svg"
# Ensure output directory exists
mkdir -p results

# STEP 1: RUN SVM PREDICTOR
echo "Running predictor..."
python3 "$PREDICTOR" "$INPUT_FASTA"
echo "Predictions written to: $OUTPUT_CSV"

# STEP 2: GENERATE FRAGMENTS USING PREDICTIONS
echo "Generating fragments..."
python3 "$FRAGMENTS" "$INPUT_FASTA" "$OUTPUT_CSV" > "$OUTPUT_FRAG"
echo "Fragments written to: $OUTPUT_FRAG"

# STEP 3: COMPUTE PHYSICOCHEMICAL PROPERTIES
echo "Computing physicochemical properties..."
python3 "$PHYCH" "$OUTPUT_FRAG" > "$OUTPUT_PHYSCH"
echo "Physicochemical properties written to: $OUTPUT_PHYSCH"

# STEP 4: PLOT SCATTERPLOT
echo "Making scatterplot..."
python3 "$SCAT" "$OUTPUT_PHYSCH"
echo "Physicochemical properties written to: $OUTPUT_SCAT"
