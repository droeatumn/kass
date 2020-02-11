# kass
main.nf assembles KIR haplotypes from PacBio HiFi reads
align.nf aligns and reports on the raw and/or assembled sequences.

<h2>Dependancies</h2>
Install Nextflow, Docker, and Git.
Create accounts in GitHub and Docker Hub.
Add 'docker.enabled = true' and 'docker.fixOwnership = true' to your Nexflow configuration (e.g., $HOME/.nextflow/config).

<h2>Assembly</h2>
<b>Input</b> <br>
The input is a directory containing one or more compressed fastq files, each representing one individual. Each file is from PacBio HiFi 99.9% consensus sequences.

<b>Output</b> <br>
Each input file has a correspoding output file (contigs.fasta)  with the assembled contigs.

<b>Running</b><br>
Use the 'base' argument to indicate the path to where kass was cloned from
GitHub, 'raw' to indicate the input directory, and 'output' to indicate the directory to put the output.

<code>    ./main.nf --base cloneDir --raw inDir --output outDir</code><br>
e.g.,
<code>    ./main.nf --base ~/git/kass --raw ~/input --output ~/output</code>

<h2>Miscellaneous</h2>
Hardware<br>
Minimum recommended hardware is 10G memory and 8 cores. More of each helps.
Run time is 1-2 hours per ID, depending on platform, genotype variation, and parallel execution.
