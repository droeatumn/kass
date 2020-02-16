#!/usr/bin/env nextflow

/*
 * Assemble KIR.
 * 
 * Requires 'docker.enabled = true' in Nexflow configuration (e.g., $HOME/.nextflow/config).
 * 
 * @author Dave Roe
 *
 */

// things that may change per run
// here are the FASTA/Q files
params.base = baseDir
base = params.base + "/"
params.raw = base + "raw/"
raw = params.raw + "/"
params.output = base + "output/"
output = params.output + "/"
params.off = 0
fqNameSuffix = "fastq.gz"          // extension on the file name (todo: expand this)
params.canuPB1 = "-pacbio-corrected"
params.canuPB2 = "-pacbio-corrected"

// things that probably won"t change per run
fqPath = raw + "*" + fqNameSuffix
markerFile = file("${baseDir}/input/markers.fasta") // gene markers
markerCapFile = file("${baseDir}/input/markers_wCap.fasta") // gene markers + capture probes
haps = base + "${baseDir}/input/HapSet23_v1.txt"
maxMem = "24g"
params.container = "droeatumn/kass:latest"

fqs = Channel.fromPath(fqPath).ifEmpty { error "cannot find any files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }

/*
 * extract
 * 
 * Extract the KIR fastq reads from potentially larger set of reads.
 * Base and output is fastq
 * Reads that don't match any kmer go to <name>_off-kir.fastq.
 * All others go to <name>_kir.fastq.
 *
 */
process extract {
  container = 'droeatumn/kass:latest'
  publishDir output, mode: 'copy', overwrite: true
  input:
    set s, file(fa) from fqs
    file(markerCapFile)
  output:
	set s, file{"*_kir.fastq"} into kirFastqs
    set s, file{"*_off-kir.fastq.gz"} into offkirFastqs optional true

  script:
    offFile="${s}_off-kir.fastq"
    offStr="out=${offFile}"
    if((params.off == null) || (params.off == 0)) {
        offStr=""
    }
	// todo: take out rname
    """
    bbduk.sh in=${fa} ${offStr} outm=${s}_kir.fastq ref=${markerCapFile} k=25 maskmiddle=f overwrite=t rename=t nzo=t rcomp=t ignorebadquality=t -Xmx${maxMem}
    find . -type f -size 0 -print0 |xargs -0 rm -f

    if [ -f ${offFile} ]; 
    then
        gzip ${offFile}
    fi
    #todo gzip ${s}_off-kir.fastq
    """
} // extract

/*
 * correct
 * 
 * Error correct the fastq reads.
 *
 */
process correct {
  container = 'droeatumn/kass:latest'
  publishDir output, mode: 'copy', overwrite: true
  input:
    set s, file(fq) from kirFastqs
  output:
	set s, file{"${s}-corrected.fasta.gz"} into correctedReads
	
    """
    lorma.sh ${fq}
    mv final.fasta ${s}-corrected.fasta
    gzip ${s}-corrected.fasta
    """
} // correct

/*
 * assemble
 * 
 * With two files per sample.
 * e.g., Kir5_P2_10_999F.fasta.gz and Kir5_P2_10_999F_gap.fasta.gz
 * 
 * 
 */
process assemble {
  container = 'droeatumn/kass:latest'
  publishDir output, mode: 'copy', overwrite: true //todo
  input:
    set s, file(fq) from correctedReads
  output:
    set s, file{"${s}*.contigs.fasta"} into assembly
	
    """
    canu -p ${s} -d ${s} rawErrorRate=0.05 correctedErrorRate=0.01 genomeSize=200k "batOptions=-dg 0.05 -db 0.05 -dr 0.05 -ca 500 -cp 50" ${params.canuPB2} ${s}-corrected.fasta.gz
    cp ${s}/${s}.contigs.fasta .
    """
} // assemble

// get the per-sample name
def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  int end = name.indexOf(fqNameSuffix)
  if ( end <= 0 ) {
    throw new Exception( "Expected file " + name + " to end in '" + fqNameSuffix + "'" );
  }
  end = end -1 // Remove the trailing '.'
  return name.substring(start, end)
} // sample

workflow.onComplete {
  println "DONE: ${ workflow.success ? 'OK' : 'FAILED' }"
}

workflow.onError {
  println "ERROR: ${workflow.errorReport.toString()}"
}

