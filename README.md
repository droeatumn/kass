# kass
main.nf assembles KIR haplotypes from PacBio HiFi reads.<br>
align.nf aligns and reports on the raw and/or assembled sequences.

<h2>Dependancies</h2>
Install Nextflow, Docker, and Git.
Create accounts in GitHub and Docker Hub.
Add 'docker.enabled = true' and 'docker.fixOwnership = true' to your Nexflow
configuration (e.g., $HOME/.nextflow/config). Make sure Docker is running
and you are logged in to Docker Hub.

<h2>Structural analysis</h2>
KPI can be used to determine the presence/absence of genes and haplotype pairs from the raw data: https://github.com/droeatumn/kpi

<h2>Assembly</h2>
<b>Input</b>
The input is a directory containing one or more compressed fastq files, each representing one individual. Each file is from PacBio HiFi 99.9% consensus sequences.<br>

<b>Output</b>
Each input file has a correspoding output file (contigs.fasta)  with the assembled contigs.

<b>Running</b>
Use the 'base' argument to indicate the path to where kass was cloned from
GitHub, 'raw' to indicate the input directory, and 'output' to indicate the directory to put the output.

<code>    ./main.nf --base cloneDir --raw inDir --output outDir</code><br>
e.g.,
<code>    ./main.nf --base ~/git/kass --raw ~/input --output ~/output</code>

The image contains an example: a cA01&tilde;tA01 homozygous individual (GenBank accession KP420442). <br>
<code>    ./main.nf --base cloneDir --raw ~/git/kass/input/example1 --output outDir</code>

<h2>Alignment</h2>
<b>Input</b>
The input is a directory containing a reference sequence in a fasta file along with one or more fasta/fastq files to be aligned to that reference.

<b>Output</b> <br>
Index files are output for the reference fasta.<br>
For each non-reference input file, a sorted bam file, its index, and the unaligned reads are output. Also, Qualimap (qualimap.pdf) and NanoPlot (NanoPlot-report.html) reports are generated for the alignment and a FastQC report (fastqc.html) is generated if the input is a fastq file.

<b>Running</b><br>
Use the 'base' argument to indicate the path to where kass was cloned from
GitHub, 'raw' to indicate the input directory, and 'output' to indicate the directory to put the output. Use 'refFasta' to indicate the name of the reference fasta file that is located in the input directory. Use 'threadNum' to optionally set maximum number of threads to use.

<code>    align.nf --base cloneDir --raw inDir --reference refFasta --output outDir --threads threadNum</code><br>
e.g.,
<code>    align.nf --raw ~/input --reference KP420442.fasta --output ~/output --threads 12</code>


<h2>Bundled references</h2>
Some references and their indexes are bundled in input/references/. 
KP420439 and KP420442 are cA01&tilde;tA01. KP420440 is cB01&tilde;tB01.

<h2>Miscellaneous</h2>
<b>Hardware</b><br>
Minimum recommended hardware is 10G memory and 8 cores. More of each helps.
Run time is 1-2 hours per ID, depending on platform, genotype variation, and parallel execution.
