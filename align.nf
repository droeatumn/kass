#! /usr/bin/env nextflow

/*
 * Aligns fasta/fastq files to a reference and reports on the alignments and raw fastq file.
 *
 * Uses bowtie for alignment.
 * Generates QualiMap and Nanoplot reports for alignments.
 * Generates FastQC reports for fastq files.
 *
 * usage : ./align.nf --threads # --reference <fasta ref file name>
 * 
 * 
 * @author Dave Roe
 * @todo add support for bwa
 */

params.home = baseDir
home = params.home + "/"
params.raw = "/opt/kass/raw"
params.output = home + "/output"
output = params.output + "/"
params.reference = ""
params.bowtie = "--end-to-end -N0"
params.threads = "8"
params.container = "droeatumn/kass:latest"
params.nocontainer = "null"

refFile = file("${params.raw}/${params.reference}")

raw = "${params.raw}/*{fastq,fq,fastq.gz,fq.gz,fasta,fa,fasta.gz,fa.gz}"
reads = Channel.fromPath(raw).ifEmpty { error "cannot find any reads matching ${raw}" }.map { path -> tuple(sample(path), path) }
alignmentReads = Channel.create()
readTap = reads.tap(alignmentReads).filter{ it[1] != refFile }
//System.err.println readTap

process align {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  // r is the input to be aligned
  input:
    set s, file(r) from readTap
    file ref from refFile
  output:
//    set s, file {"*_sorted.bam*"} into bamOutput
    tuple val(s), file("*_sorted.bam"), file("*_sorted.bam.bai"), file(r) into bamOutput
    set s, file {"*_unaligned.txt.gz"} into unaligned
    set s, file {"*.fasta.*" } into refs
    set s, file {"bowtie_indexes" } into btIndexes
    
  script:
    // todo: modularize this chunk
    def refName = ref.name.replaceFirst(".fasta", "")
    def bamName = r.name.replaceFirst(".fastq.gz", "").replaceFirst(".fasta.gz", "").replaceFirst(".fastq", "").replaceFirst(".fasta", "")
    def outName = bamName + "_" + refName
    def fastaFlag = ""  // tell bowtie it is a fasta or fastq
    if((r.name.endsWith("fasta") || r.name.endsWith("fa")) ||
       r.name.endsWith("fasta.gz") || r.name.endsWith("fa.gz")) {
        fastaFlag = "-f"
    }
    
    """
    # indexing
    # todo: split this; only needs to happen once
    echo ${s} ${r} ${params.reference} ${refName} $fastaFlag
    samtools faidx ${params.reference}
    bwa index ${params.reference}
    mkdir -p bowtie_indexes
    cd bowtie_indexes
    bowtie2-build ../${params.reference} ${refName}_index
    cd ..

    # alignment  
    # bowtie2 ${params.bowtie} --threads ${params.threads} -x bowtie_indexes/${refName}_index ${fastaFlag} ${r} -S ${outName}.sam --un ${outName}_unaligned.txt 2> ${outName}_err.txt
    bwa mem -xpacbio -t${params.threads} ${ref} ${r} > ${outName}.sam 2> ${outName}_err.txt
    samtools view -b -S ${outName}.sam > ${outName}.bam
    samtools sort ${outName}.bam -o ${outName}_sorted.bam
    samtools index ${outName}_sorted.bam
    rm ${outName}.[bs]am
    samtools view -b -f 4 ${outName}_sorted.bam > ${outName}_unaligned.txt
    gzip -f ${outName}_unaligned.txt

    """
} // align

process report {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    errorStrategy 'ignore'
    publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  // r is the input to be aligned
  input:
    set s, file(bamFile), file(bamIndexFile), file(r) from bamOutput
    file ref from refFile
  output:
    set s, file {"*_reports" } into reports

  script:
    // todo: modularize this chunk
    def refName = ref.name.replaceFirst(".fasta", "")
    def bamName = r.name.replaceFirst(".fastq.gz", "").replaceFirst(".fasta.gz", "").replaceFirst(".fastq", "").replaceFirst(".fasta", "")
    def outName = bamName + "_" + refName
    def fastaFlag = ""  // tell bowtie it is a fasta or fastq
    if((r.name.endsWith("fasta") || r.name.endsWith("fa")) ||
       r.name.endsWith("fasta.gz") || r.name.endsWith("fa.gz")) {
        fastaFlag = "-f"
    }

    """
    mkdir -p ${outName}_reports
    qualimap bamqc -bam ${outName}_sorted.bam -gd HUMAN -outdir ${outName}_reports -outformat PDF
    mv ${outName}_reports/report.pdf ${outName}_reports/${outName}_qualimap.pdf
    mv ${outName}_reports/genome_results.txt ${outName}_reports/${outName}_qualimap_genome_results.txt
    NanoPlot -t ${params.threads} --bam ${outName}_sorted.bam -o ${outName}_reports -p ${outName} -f pdf --N50
    if ["$fastaFlag" == ""]
    then
        fastqc -t ${params.threads} -o ${outName}_reports -f fastq ${r}
    fi
    """
} // report

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf('.'))
}
