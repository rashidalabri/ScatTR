#!/usr/bin/env bash

# Make sure we're in the right directory
scriptDir=$(cd $(dirname $0) && pwd -P)
cd $scriptDir

# Set up output directory and prefix
OUTPUT_DIR=output
OUTPUT_PREFIX=$OUTPUT_DIR/sample

mkdir -p $OUTPUT_DIR

if command -v cargo &> /dev/null; then
    RUN_CMD="cargo run --quiet --release --"
elif command -v scattr &> /dev/null; then
    RUN_CMD="scattr"
else
    echo "Neither 'cargo' nor 'scattr' found in PATH."
    exit 1
fi

# Run ScatTR end-to-end

# 1. Extract the bag of reads 
# Output: output/sample.bag.bam
$RUN_CMD $OUTPUT_PREFIX extract \
    input/sample.bam \
    input/catalog.tsv 

# 2. Extract insert and depth distributions
# Output: output/sample.stats.json
$RUN_CMD $OUTPUT_PREFIX stats \
    input/sample.bam \
    --include-contigs test_genome

# 3. Define the optimization problems
# Output: output/sample.defs.json
$RUN_CMD $OUTPUT_PREFIX define \
    input/catalog.tsv \
    input/reference.fa

# 4. Genotype the sample
# Output: output/sample.genotypes.json
$RUN_CMD $OUTPUT_PREFIX genotype