#! /usr/bin/env nextflow

/*
 * Annotates strings in fasta files.
 * 
 * @author Dave Roe
 * @todo handle gzipped input
 * @todo standardize "KIR" in the output
 * @todo remove the bar in the ft output
 */

params.home = baseDir
home = params.home + "/"
params.raw = home + '/raw/'
//params.output = "/Users/daver/git/kass/output"
params.output = home + "/output"
output = params.output + "/"
params.bowtie = "-a --end-to-end --rdg 3,3 --rfg 3,3"
params.threads = "28"
params.maxMem = "200G"
params.container = "droeatumn/kass:latest"
params.nocontainer = "null"
params.sbt = null

refAlleleDir = file("${home}/input/references/genes")
jointAlleleDir = file("${home}/input/references/joints")
featuresFile = file("${home}/input/features.txt") // markup features
markerFile = file("${home}/input/cap.fasta") // capture markers
alignProbesFile = file("${home}/src/alignment2ProbePairs.groovy")
annotateFile = file("${home}/src/annotateMarkup.groovy")
if(params.sbt == null) { 
    sbtFile = home + "/input/SubmissionTemplate.sbt"
} else { 
    sbtFile = file(params.raw + "/" + params.sbt)
}

raw = "${params.raw}/*{fasta,fa,fasta.gz,fa.gz}"
reads = Channel.fromPath(raw).ifEmpty { error "cannot find any reads matching ${raw}" }.map { path -> tuple(sample(path), path) }
annotateReads = Channel.create()
readTap = reads.tap(annotateReads).filter{ it[1] != params.sbt }

/* 
 * Use the capture probes to orient the input sequences.
 */
process orient {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    tag { s }

  input:
    set s, file(r) from readTap
    path(markerFile)
  output:
    tuple s, file{"${s}*-orient.fasta"} into orientedFasta
  script:
    def rootName = r.name.replaceFirst(".gz", "")
    def gzFlag = ""  // tell bowtie it is a fasta or fastq
    if(r.name.endsWith(".gz")) {
        gzFlag = "1"
    }
    """
    if [ "$gzFlag" == "1" ]; then
        gunzip -f ${r}
    fi
    orient.groovy -i ${rootName} -p ${markerFile} -o tmp.fasta
    reformat.sh -Xmx${params.maxMem} in=tmp.fasta out=${s}-orient.fasta overwrite=true
    """
} // orient

/* 
 * Use the capture probes to annotate.
 * Output is *_annotation.txt
 */
process structure {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    tag { s }

  input:
    set s, file(r) from orientedFasta
    path(markerFile)
    path(alignProbesFile)
    path(annotateFile)
    path(featuresFile)
  output:
    tuple s, file{"${s}*_markup.txt"} into markup
    tuple s, path{"${s}*_annotation.txt"} into annotation
    tuple s, path{r}, path{"${s}*_features.fasta"} into features mode flatten
    
  script:
    // todo: modularize this chunk
    def bamName = r.name.replaceFirst(".fastq.gz", "").replaceFirst(".fasta.gz", "").replaceFirst(".fastq", "").replaceFirst(".fasta", "").replaceFirst("-orient", "")
    def outName = bamName
    def fastaFlag = ""  // tell bowtie it is a fasta or fastq
    if((r.name.endsWith("fasta") || r.name.endsWith("fa")) ||
       r.name.endsWith("fasta.gz") || r.name.endsWith("fa.gz")) {
        fastaFlag = "-f"
    }
    
    """
    # indexing
    # todo: split this; only needs to happen once
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
    echo "done"
    deep.pl replace '2DL4~3DL2' '2DL4~3DL1L2' '*features.fasta' || TRUE
    deep.pl replace '2DL4~3DL2' '2DL4~3DL1S1' '*_annotation.txt' || TRUE
    cut -f2 ${bamName}_annotation.txt | sort | uniq -c > ${bamName}_annotation_strings.txt

    """
} // structure

process interpret {
    if(params.nocontainer == "null") {
        container = params.container
    }
    //publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  input:
    set s, file(inContig), file(r) from features
    path(refAlleleDir)
  output:
    path{"${s}*.tbl.txt"} optional true into tables
    
  script:
    def nameNoExt = r.name.replaceAll("_features.fasta", "")
    def intervening = nameNoExt.contains("3DP1-2DL4") // skip this one
    def i = nameNoExt.lastIndexOf('_')
    def locus = nameNoExt[i+1..-1]

    """
    if [ "${locus}" != "3DP1-2DL4" ]; then
        gfi.groovy -i ${refAlleleDir} -f ${r} -g KIR${locus} -j ${jointAlleleDir} -o . 2> gfi_${s}_${locus}_err.txt
    fi
    """

} // 

allFTs = tables.collectFile(name: 'tbl.txt', newLine: true)

/*
process publishGL {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    input:
        path(g) from gl
    output:
        path(g)
    script:
    """
    """
} // publishGL
*/

/*
 * Combine the per-feature feature tables into the per contig annotation.
 * @todo linux64.table2asn_GFF -augustus-fix -f ${desc}.gff -i ${inContig} -outdir . -genbank -verbose -euk -V b -Z  -t ${sbtF} -j "[organism=Homo sapiens]"

 */
process combine {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    input:
        path(ft) from allFTs
    output:
        path("*.ft.txt")
    script:
    """
    echo "${ft}"
    combineFT.groovy -i tbl.txt 2> combineFT_err.txt
    """
} // combine

// combineGL (todo)

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf('.'))

}
