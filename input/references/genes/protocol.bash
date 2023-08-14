#!/bin/bash
# To create KIR2DP1_prot.fasta and KIR3DP1_prot.fasta, load the _gen.msf
# files into a protein translation tools (e.g.,
# https://www.bioinformatics.org/sms2/translate.html). 
#

FASTA_DIR=$HOME/git/IPDKIR/fasta/
DEEP=$HOME/git/kass/bin/deep.pl

#
# nuc
cp $FASTA_DIR/KIR[23]*_nuc.fasta .
$DEEP replace '*' '_' '*nuc.fasta'
$DEEP replace 'KIR:KIR\d+ ' '' '*nuc.fasta' --literal=0
$DEEP replace ' .*' '' '*nuc.fasta' --literal=0
cat KIR3DL1_nuc.fasta KIR3DS1_nuc.fasta > KIR3DL1S1_nuc.fasta
cat KIR2DL2_nuc.fasta KIR2DL3_nuc.fasta > KIR2DL2L3_nuc.fasta
cat KIR2DL1_nuc.fasta KIR2DS1_nuc.fasta KIR2DS2_nuc.fasta > KIR2DL1S1S2_nuc.fasta
cat KIR2DS3_nuc.fasta KIR2DS4_nuc.fasta KIR2DS5_nuc.fasta > KIR2DS3S4S5_nuc.fasta
cat KIR2DL2L3_nuc.fasta KIR2DS3S4S5_nuc.fasta > KIR2DL2L3S3S4S5_nuc.fasta

# protein
cp $FASTA_DIR/KIR[23]*_prot.fasta .
$DEEP replace '*' '_' '*prot.fasta'
$DEEP replace 'KIR:KIR\d+ ' '' '*prot.fasta' --literal=0
$DEEP replace ' .*' '' '*prot.fasta' --literal=0
cat KIR3DL1_prot.fasta KIR3DS1_prot.fasta > KIR3DL1S1_prot.fasta
cat KIR2DL2_prot.fasta KIR2DL3_prot.fasta > KIR2DL2L3_prot.fasta
cat KIR2DL1_prot.fasta KIR2DS1_prot.fasta KIR2DS2_prot.fasta > KIR2DL1S1S2_prot.fasta
cat KIR2DS3_prot.fasta KIR2DS4_prot.fasta KIR2DS5_prot.fasta > KIR2DS3S4S5_prot.fasta
cat KIR2DL2L3_prot.fasta KIR2DS3S4S5_prot.fasta > KIR2DL2L3S3S4S5_prot.fasta

# gene
cp $FASTA_DIR/KIR[23]*_gen.fasta .
$DEEP replace '*' '_' '*gen.fasta'
$DEEP replace 'KIR:KIR\d+ ' '' '*gen.fasta' --literal=0
$DEEP replace ' .*' '' '*gen.fasta' --literal=0
cat KIR3DL1_gen.fasta KIR3DS1_gen.fasta > KIR3DL1S1_gen.fasta
cat KIR2DL2_gen.fasta KIR2DL3_gen.fasta > KIR2DL2L3_gen.fasta
cat KIR2DL1_gen.fasta KIR2DS1_gen.fasta KIR2DS2_gen.fasta > KIR2DL1S1S2_gen.fasta
cat KIR2DS3_gen.fasta KIR2DS4_gen.fasta KIR2DS5_gen.fasta > KIR2DS3S4S5_gen.fasta
cat KIR2DL2L3_gen.fasta KIR2DS3S4S5_gen.fasta > KIR2DL2L3S3S4S5_gen.fasta

rm KIR2DL1_*.fasta KIR2DL2L3_*.fasta KIR2DL2_*.fasta KIR2DL3_*.fasta KIR2DS1_*.fasta KIR2DS2_*.fasta KIR2DS3_*.fasta KIR2DS5_*.fasta KIR3DL1_*.fasta KIR3DS1_*.fasta
