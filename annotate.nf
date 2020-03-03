#! /usr/bin/env nextflow

/*
 * Annotates strings in fasta files.
 * 
 * @author Dave Roe
 * 
 */

params.home = baseDir
home = params.home + "/"
params.raw = "/opt/kass/raw"
//params.output = "/Users/daver/git/kass/output"
params.output = home + "/output"
output = params.output + "/"
params.reference = ""
params.bowtie = "-a --end-to-end --rdg 3,3 --rfg 3,3"
params.threads = "8"
params.container = "droeatumn/kass:latest"
params.nocontainer = "null"

refFile = file("${params.raw}/${params.reference}")
featuresFile = file("${home}/input/features.txt") // markup features
markerFile = file("${home}/input/markers.fasta") // gene markers
alignProbesFile = file("${home}/src/alignment2ProbePairs.groovy")
annotateFile = file("${home}/src/annotateMarkup.groovy")

raw = "${params.raw}/*{fasta,fa,fasta.gz,fa.gz}"
reads = Channel.fromPath(raw).ifEmpty { error "cannot find any reads matching ${raw}" }.map { path -> tuple(sample(path), path) }
alignmentReads = Channel.create()
readTap = reads.tap(alignmentReads).filter{ it[1] != refFile }
//System.err.println readTap

process annotateStructure {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  // r is the input to be aligned
  input:
    set s, file(r) from readTap
    file ref from refFile
    path(markerFile)
    path(alignProbesFile)
    path(annotateFile)
    path(featuresFile)
  output:
//    set s, file {"*_sorted.bam*"} into bamOutput
//    tuple val(s), file("*_sorted.bam"), file("*_sorted.bam.bai"), file(r) into bamOutput
    tuple s, file{"${s}*_markup.txt"} into markup
    tuple s, path{"${s}*_annotation.txt"} into annotation
    tuple s, path{"${s}*_annotation_strings.txt"} into annotationStrings
    
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
    echo ${s} ${r} ${refName} $fastaFlag
    echo "indexing ${r}..."
    samtools faidx ${r}
    #bwa index ${r}
    mkdir -p bowtie_indexes
    cd bowtie_indexes
    bowtie2-build --threads ${params.threads} ../${r} ${s}_index
    cd ..

    # alignment  
    echo "aligning ${markerFile} to ${r}..."
    bowtie2 ${params.bowtie} --threads ${params.threads} -x bowtie_indexes/${s}_index ${fastaFlag} ${markerFile} -S ${outName}.sam 2> ${outName}_err.txt
    samtools view -b -S ${outName}.sam > ${outName}.bam
    #samtools sort ${outName}.bam -o ${outName}_sorted.bam
    #samtools index ${outName}_sorted.bam
    samtools sort ${outName}.bam -O sam -o ${outName}_sorted.sam
    rm ${outName}.[bs]am

    # markup the alignment of the probe pairs
    mkdir -p annotation
    echo "marking ${bamName}..."
    ./${alignProbesFile} -d . -m 1 -o ${bamName}_markup.txt 2> ${bamName}_markup_err.txt
    tail ${bamName}_markup_err.txt
    
    echo "annotating ${s}..."
    # annotate the markup with the genes
    ./${annotateFile} -i ${featuresFile} -f ${r} -m ${bamName}_markup.txt -o . 2> ${bamName}_annotation_err.txt
    cut -f2 ${bamName}_annotation.txt | sort | uniq -c > ${bamName}_annotation_strings.txt

    """
} // annotateStructure

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf('.'))
}
