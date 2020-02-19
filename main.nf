#!/usr/bin/env nextflow

/*
 * Assemble KIR.
 * 
 * Requires 'docker.enabled = true' in Nexflow configuration (e.g., $HOME/.nextflow/config).
 * 
 * @author Dave Roe
 * @todo --resume isn't working
 */

// things that may change per run
// here are the FASTA/Q files
params.base = baseDir
base = params.base + "/"
params.raw = base + "/raw/"
raw = params.raw + "/"
params.output = base + "output/"
output = params.output + "/"
params.off = 0
fqNameSuffix = "fastq.gz"          // extension on the file name (todo: expand this)
params.canuPB1 = "-pacbio-corrected"
params.canuPB2 = "-pacbio-corrected"

// things that probably won"t change per run
fqPath = raw + "/*/" + fqNameSuffix
markerFile = file("${baseDir}/input/markers.fasta") // gene markers
markerCapFile = file("${baseDir}/input/markers_wCap.fasta") // gene markers + capture probes
featuresFile = file("${baseDir}/input/features.txt") // markup features
haps = base + "${baseDir}/input/HapSet23_v1.txt"
alignProbesFile = file("${baseDir}/src/alignment2ProbePairs.groovy")
annotateFile = file("${baseDir}/src/annotateMarkup.groovy")
maxMem = "24g"
params.container = "droeatumn/kass:latest"

fqs = Channel.fromPath(fqPath).ifEmpty { error "cannot find any files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }

/*
 * extract
 * 
 * Extract the KIR fastq reads from potentially larger tuple of reads.
 * Base and output is fastq
 * Reads that don't match any kmer go to <name>_off-kir.fastq.
 * All others go to <name>_kir.fastq.
 *
 */
process extract {
  container = 'droeatumn/kass:latest'
  publishDir output, mode: 'copy', overwrite: true
  input:
    tuple s, path(fa) from fqs
    path(markerCapFile)
  output:
	tuple s, file{"*_kir.fastq"} into kirFastqs
    tuple s, file{"*_off-kir.fastq.gz"} into offkirFastqs optional true

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
    tuple s, path(fq) from kirFastqs
  output:
	tuple s, file{"${s}-corrected.fasta.gz"} into correctedReads
	
    """
    lorma.sh ${fq}
    mv final.fasta ${s}-corrected.fasta
    gzip ${s}-corrected.fasta
    """
} // correct

/*
 * assemble
 * 
 */
process assemble {
  container = 'droeatumn/kass:latest'
  publishDir output, mode: 'copy', overwrite: true //todo
  input:
    tuple s, path(fq) from correctedReads
  output:
    tuple s, file{"${s}*.contigs.fasta"} into assembly
	
    """
    canu -p ${s} -d ${s} rawErrorRate=0.05 correctedErrorRate=0.01 genomeSize=200k "batOptions=-dg 0.05 -db 0.05 -dr 0.05 -ca 500 -cp 50" ${params.canuPB2} ${s}-corrected.fasta.gz
    cp ${s}/${s}.contigs.fasta .
    """
} // assemble

process annotateStructure {
  container = 'droeatumn/kass:latest'
  publishDir output, mode: 'copy', overwrite: true

  input:
    tuple s, path(contigs) from assembly
    path(markerFile)
  output:
    tuple s, file{"${s}*_markup.txt"} into markup
    tuple s, path{"${s}*_annotation.txt"} into annotation
	
    """
    # make indexes
    samtools faidx ${contigs}
    mkdir -p bowtie_indexes
    cd bowtie_indexes
    bowtie2-build ../${contigs} ${s}
    cd ..

    # align
    nice bowtie2 -a --end-to-end --rdg 3,3 --rfg 3,3 -p8 -x bowtie_indexes/${s} -f ${markerFile} -S ${s}.sam 2> ${s}_err.txt
    cat ${s}_err.txt
    samtools view -b -S ${s}.sam > ${s}.bam
    samtools sort ${s}.bam -o ${s}_sorted.bam
    samtools index ${s}_sorted.bam
    samtools sort ${s}.bam -O sam -o ${s}_sorted.sam
    rm ${s}.[bs]am

    # markup the alignment of the probe pairs
    mkdir -p annotation
    ${alignProbesFile} -d . -m 1 -o ${s}_markup.txt 2> ${s}_markup_err.txt
    tail ${s}_markup_err.txt
    
    # annotate the markup with the genes
    ${annotateFile} -i ${featuresFile} -f ${contigs} -m ${s}_markup.txt -o . 2> ${s}_annotation_err.txt
    cut -f2 ${s}_annotation.txt | sort | uniq -c > ${s}_annotation_strings.txt
    """
} // annotate

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

