# kass
main.nf assembles KIR haplotypes from PacBio HiFi reads
align.nf aligns and reports on the raw and/or assembled sequences.

<h2>Dependancies</h2>
Install Nextflow, Docker, and Git.
Create accounts in GitHub and Docker Hub.
Add 'docker.enabled = true' and 'docker.fixOwnership = true' to your Nexflow configuration (e.g., $HOME/.nextflow/config). Make sure Docker is running.

<h2>Assembly</h2>
<b>Input</b>
The input is a directory containing one or more compressed fastq files, each representing one individual. Each file is from PacBio HiFi 99.9% consensus sequences.

<b>Output</b>
Each input file has a correspoding output file (contigs.fasta)  with the assembled contigs.

<b>Running</b>
Use the 'base' argument to indicate the path to where kass was cloned from
GitHub, 'raw' to indicate the input directory, and 'output' to indicate the directory to put the output.

<code>    ./main.nf --base cloneDir --raw inDir --output outDir</code><br>
e.g.,
<code>    ./main.nf --base ~/git/kass --raw ~/input --output ~/output</code>

The image contains an example: a cA01&tilde;tA01 homozygous individual (GenBank accession KP420442). <br>
<code>    ./main.nf --base ~/git/kass --raw ~/git/kass/input/example1 --output ~/output</code>

<h2>Bundled references</h2>
Some references and their indexes are bundled in input/references/. 
KP420439 and KP420442 are cA01&tilde;tA01. KP420440 is cB01&tilde;tB01.

<h2>Miscellaneous</h2>
Hardware<br>
Minimum recommended hardware is 10G memory and 8 cores. More of each helps.
Run time is 1-2 hours per ID, depending on platform, genotype variation, and parallel execution.
