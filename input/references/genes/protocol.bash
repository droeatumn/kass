#!/bin/bash
# To create KIR2DP1_prot.fasta and KIR3DP1_prot.fasta, load the _gen.msf
# files into Aliview. Then File -> Save as Translated alignment (Amino Acid).
#
# Align KIR2DP1_prot.fasta and KIR3DP1_prot.fasta with MAFFT (v7).
# https://mafft.cbrc.jp/alignment/server/index.html
# Align the combined loci (e.g., KIR3DL1S1_prot_msa.fasta) by aligning the combined fasta (e.g., KIR3DL1S1_prot.fasta).
# 
# To convert _prot.msf files from msf format to fasta (e.g.
# KIR2DL1_prot_msa.fasa), use Aliview or the seqret web server
# https://www.ebi.ac.uk/Tools/sfc/emboss_seqret/
#
# Run this script, then blat_ipd.bash, hints_ipd.bash, proteinprofile_ipd.bash.

#
# nuc
cp ../fasta/KIR[23]*_nuc.fasta .
deep.pl replace '*' '_' '*nuc.fasta'
deep.pl replace 'IPD:KIR\d+ ' '' '*nuc.fasta' --literal=0
deep.pl replace ' .*' '' '*nuc.fasta' --literal=0
cat KIR3DL1_nuc.fasta KIR3DS1_nuc.fasta > KIR3DL1S1_nuc.fasta
cat KIR2DL2_nuc.fasta KIR2DL3_nuc.fasta > KIR2DL2L3_nuc.fasta
cat KIR2DL1_nuc.fasta KIR2DS1_nuc.fasta KIR2DS2_nuc.fasta > KIR2DL1S1S2_nuc.fasta
cat KIR2DS3_nuc.fasta KIR2DS4_nuc.fasta KIR2DS5_nuc.fasta > KIR2DS3S4S5_nuc.fasta
cat KIR2DL2L3_nuc.fasta KIR2DS3S4S5_nuc.fasta > KIR2DL2L3S3S4S5_nuc.fasta

# protein
cp ../fasta/KIR[23]*_prot.fasta .
deep.pl replace '*' '_' '*prot.fasta'
deep.pl replace 'IPD:KIR\d+ ' '' '*prot.fasta' --literal=0
deep.pl replace ' .*' '' '*prot.fasta' --literal=0
cat KIR3DL1_prot.fasta KIR3DS1_prot.fasta > KIR3DL1S1_prot.fasta
cat KIR2DL2_prot.fasta KIR2DL3_prot.fasta > KIR2DL2L3_prot.fasta
cat KIR2DL1_prot.fasta KIR2DS1_prot.fasta KIR2DS2_prot.fasta > KIR2DL1S1S2_prot.fasta
cat KIR2DS3_prot.fasta KIR2DS4_prot.fasta KIR2DS5_prot.fasta > KIR2DS3S4S5_prot.fasta
cat KIR2DL2L3_prot.fasta KIR2DS3S4S5_prot.fasta > KIR2DL2L3S3S4S5_prot.fasta

# gene
cp ../fasta/KIR[23]*_gen.fasta .
deep.pl replace '*' '_' '*gen.fasta'
deep.pl replace 'KIR:KIR\d+ ' '' '*gen.fasta' --literal=0
deep.pl replace ' .*' '' '*gen.fasta' --literal=0
cat KIR3DL1_gen.fasta KIR3DS1_gen.fasta > KIR3DL1S1_gen.fasta
cat KIR2DL2_gen.fasta KIR2DL3_gen.fasta > KIR2DL2L3_gen.fasta
cat KIR2DL1_gen.fasta KIR2DS1_gen.fasta KIR2DS2_gen.fasta > KIR2DL1S1S2_gen.fasta
cat KIR2DS3_gen.fasta KIR2DS4_gen.fasta KIR2DS5_gen.fasta > KIR2DS3S4S5_gen.fasta
cat KIR2DL2L3_gen.fasta KIR2DS3S4S5_gen.fasta > KIR2DL2L3S3S4S5_gen.fasta
