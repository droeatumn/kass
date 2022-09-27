# kass
main.nf assembles KIR haplotypes from PacBio HiFi reads.<br>
annotate.nf annotates the structure of assembled contigs<br>
align.nf aligns and reports on the raw or assembled sequences.

<h2>Dependencies</h2>
Install Java, Groovy, Nextflow, Docker, and Git.
Create accounts in GitHub and Docker Hub.
Add 'docker.enabled = true' and 'docker.fixOwnership = true' to your Nexflow
configuration (e.g., $HOME/.nextflow/config). Make sure Docker is running
and you are logged in to Docker Hub.
Use the --nocontainer option to run without any container (natively).

<h2>Structural analysis</h2>
KPI can be used to determine the presence/absence of genes and haplotype pairs from the raw data: https://github.com/droeatumn/kpi

<h2>Assembly</h2>
<b>Input</b>
The input is a directory containing one or more compressed fastq files, each representing one individual. Each file is from PacBio HiFi consensus sequences (preferably 99.9%).<br>
<br>
<b>Output</b><br>
Each input file has a correspoding output file (*.contigs.fasta)  with the assembled contigs. Each contig is annotated with gene content and order in the file with the suffix 'annoation.txt'.<br>
<br>
<b>Running</b>
Use the parameter 'raw' to indicate the input directory, and 'output' to indicate the directory to put the output.  The defaults are 'raw' and 'output' under the location where kass was pulled. Optionally use --threads to optionally set maximum number of threads to use (default 8). To output the off-kir reads, use --off.<br>

<code>    ./main.nf --raw inDir --output outDir</code><br>
e.g.,
<code>    ./main.nf  --raw ~/input --output ~/output</code>

The image contains an example: simulated reads from a single cA01&tilde;tA01 haplotype (GenBank accession KP420442). <br>
<code>    ./main.nf --raw ~/git/kass/input/example1 --output outDir</code>

<h2>Annotation</h2>
<b>Input</b><br>
The input is a folder containing fasta files (usually contigs) to be annotated. Each file may contain more than one sequence. The files cannot be gzipped.<br>
<br>
<b>Output</b><br>
For each gene, feature tables (ft.txt) and a features sequences (features.fasta) are created to annotate the alleles. The feature table contains the gene and intron/exon locations for each contig in GenBank update format. The features.fasta files extract the alleles from the contigs and group them into their genes.<br><br>
For more information, an annotation file (annotation.txt) will be created to annotate the gene locations in each contig. Ideally the annotation only contains gene names, as opposed to markup. Markup means that part of the assembly didn't match the markup for any full gene. The markup.txt files have more information about the locations of the probes in the contigs and the key to converting the probe pairs into motif characters. If the annotation.txt files have motif markup (characters) instead of gene names, then either the assembly is wrong in that area or the regex expressions that map the motifs to the genes needs to be altered. This is defined in input/features.txt.<br>
<br>
<b>Running</b><br>
Use the 'raw' parameter to indicate the input directory, and 'output' to indicate the directory to put the output. Use 'refFasta' to indicate the name of the reference fasta file that is located in the input directory. Use 'threadNum' to optionally set maximum number of threads to use (default 8).

<code>    annotate.nf --raw inDir --output outDir --threads threadNum</code><br>
e.g.,
<code>    annotate.nf --raw ~/input --output ~/output --threads 12</code>

<h2>Alignment</h2>
<b>Input</b><br>
The input is a directory containing a reference sequence in a fasta file along with one or more fasta/fastq files to be aligned to that reference.<br>
<br>
<b>Output</b><br>
Index files are output for the reference fasta.<br>
For each non-reference input file, a sorted bam file, its index, and the unaligned reads are output. Also, Qualimap (qualimap.pdf), NanoPlot (NanoPlot-report.html), and Quast (quast/report.html) reports are generated for the alignment and a FastQC report (fastqc.html) is generated if the input is a fastq file. Also for fastq files, add <code>--bwa="-xpacbio"</code> to the command.<br>
Use caution interpreting these reports. For example, Quast doesn't just report on the alignment. It reports misassembles that are really patterns of homology. This leads to erroneously high rates of mismatches. Use qualimap for the error rate.<br>
<br>
<b>Running</b><br>
Use the 'raw' parameter to indicate the input directory, and 'output' to indicate the directory to put the output. Use 'refFasta' to indicate the name of the reference fasta file that is located in the input directory. Use 'threadNum' to optionally set maximum number of threads to use (default 8).<br>

<code>    align.nf --raw inDir --reference refFasta --output outDir --threads threadNum</code><br>
e.g.,
<code>    align.nf --raw ~/input --reference KP420442.fasta --output ~/output --threads 12</code>

<h2>Bundled references</h2>
Some references and their indexes are bundled in input/references/. 
KP420439 and KP420442 are cA01&tilde;tA01. KP420440 is cB01&tilde;tB01. They each have a bed file that documents the locations of the genes.

<h2>Citation</h2>
Roe D. Efficient Sequencing, Assembly, and Annotation of Human KIR Haplotypes. Frontiers in Immunology (2020) 11:11 (https://doi.org/10.3389/fimmu.2020.582927)
<br>
  
<h2>Miscellaneous</h2>
<b>Hardware</b><br>
Minimum recommended hardware for targeted sequencing is 30G memory and 8 CPU
cores.
More of each helps, especially with WGS.
Run time is 1-2 hours per ID, depending on platform, genotype variation, and parallel execution.
<br><br>
<b>Allele annotation</b><br>
Full-gene alleles are annotated with respect to the IPD database version 2.10.0 (http://www.ebi.ac.uk/ipd/kir).<br>
Robinson J, Halliwell JA, Hayhurst JH, Flicek P, Parham P, Marsh SGE: The IPD and IPD-IMGT/HLA Database: allele variant databases Nucleic Acids Research (2015), 43:D423-431<br>
Robinson J, Malik A, Parham P, Bodmer JG, Marsh SGE: IMGT/HLA - a sequence database for the human major histocompatibility complex Tissue Antigens (2000), 55:280-287
