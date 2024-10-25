#!/usr/bin/env bash

scriptDir=$(cd $(dirname $0) && pwd -P)
cd $scriptDir

# Motif sequence and length
MOTIF="CAGATA"
MOTIF_LEN=${#MOTIF}

# Reference copy number, expansion copy number, and flanking sequence length
REF_CN=3
EXPANSION_CN=100
FLANK=100000

# Calculate reference start and end
REF_START=$FLANK
REF_END=$((REF_START + (REF_CN * MOTIF_LEN)))

INPUT_DIR=../input
TMP_DIR=.temp_dir

echo "Motif: $MOTIF"
echo "Motif length: $MOTIF_LEN"
echo "Reference copy number: $REF_CN"
echo "Expansion copy number: $EXPANSION_CN"
echo "Flanking sequence length: $FLANK"
echo "Reference start: $REF_START"
echo "Reference end: $REF_END"
echo "Input directory: $INPUT_DIR"
echo  "Temp directory: $TMP_DIR"

mkdir -p $INPUT_DIR
mkdir -p $TMP_DIR

# Create reference and expanded genome references
python make_genome.py $MOTIF $REF_CN $TMP_DIR/reference.fa -f $FLANK
python make_genome.py $MOTIF $EXPANSION_CN $TMP_DIR/expanded.fa -f $FLANK

# Index reference genome
bwa index $TMP_DIR/reference.fa
samtools faidx $TMP_DIR/reference.fa

# Simulate reads for each haplotype (normal and expanded)
art_illumina --id hap1 -ss HSXn -sam -p --noALN -i $TMP_DIR/reference.fa -rs 42 -o $TMP_DIR/reads_hap1 -f 15 -l 150 -m 450 -s 50
art_illumina --id hap2 -ss HSXn -sam -p --noALN -i $TMP_DIR/expanded.fa -rs 42 -o $TMP_DIR/reads_hap2 -f 15 -l 150 -m 450 -s 50

# Combine reads of both haplotypes
cat $TMP_DIR/reads_hap11.fq $TMP_DIR/reads_hap21.fq > $TMP_DIR/reads1.fq
cat $TMP_DIR/reads_hap12.fq $TMP_DIR/reads_hap22.fq > $TMP_DIR/reads2.fq

# Align reads
bwa mem $TMP_DIR/reference.fa $TMP_DIR/reads1.fq $TMP_DIR/reads2.fq > $TMP_DIR/sample.sam

samtools sort -o $TMP_DIR/sample.sorted.bam $TMP_DIR/sample.sam
samtools index $TMP_DIR/sample.sorted.bam

# Create a catalog.tsv file
echo -e "id\tcontig\tstart\tend\tmotif" > $TMP_DIR/catalog.tsv
echo -e "test_locus\ttest_genome\t$REF_START\t$REF_END\t$MOTIF" >> $TMP_DIR/catalog.tsv

# Move files to input directory
mkdir -p $INPUT_DIR
mv $TMP_DIR/sample.sorted.bam $INPUT_DIR/sample.bam
mv $TMP_DIR/sample.sorted.bam.bai $INPUT_DIR/sample.bam.bai
mv $TMP_DIR/reference.fa $INPUT_DIR
mv $TMP_DIR/reference.fa.fai $INPUT_DIR
mv $TMP_DIR/catalog.tsv $INPUT_DIR

rm -r $TMP_DIR

echo "Done!"