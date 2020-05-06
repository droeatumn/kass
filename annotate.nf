#! /usr/bin/env nextflow

/*
 * Annotates strings in fasta files.
 * 
 * @author Dave Roe
 * @todo: this fails when the direction is the wrong way   *
 *    maybe add a step that standardizes based on the probes?
 *    also have to break the contigs into multiple files
 * @todo handle gzipped input
 * @todo remove full paths to augustus bin and scripts
 * @todo rename gl.txt
 * @todo standardize "KIR" in the output
 */

params.home = baseDir
home = params.home + "/"
params.raw = home + '/raw/'
//params.output = "/Users/daver/git/kass/output"
params.output = home + "/output"
output = params.output + "/"
params.bowtie = "-a --end-to-end --rdg 3,3 --rfg 3,3"
params.threads = "28"
params.container = "droeatumn/kass:latest"
params.nocontainer = "null"
params.sbt = home + "/input/SubmissionTemplate.sbt"

refAlleleDir = file("${home}/input/references/genes")
featuresFile = file("${home}/input/features.txt") // markup features
markerFile = file("${home}/input/cap.fasta") // capture markers
alignProbesFile = file("${home}/src/alignment2ProbePairs.groovy")
annotateFile = file("${home}/src/annotateMarkup.groovy")

raw = "${params.raw}/*{fasta,fa,fasta.gz,fa.gz}"
reads = Channel.fromPath(raw).ifEmpty { error "cannot find any reads matching ${raw}" }.map { path -> tuple(sample(path), path) }
annotateReads = Channel.create()
readTap = reads.tap(annotateReads).filter{ it[1] != params.sbt }

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
    set s, file(r) from readTap
    path(markerFile)
    path(alignProbesFile)
    path(annotateFile)
    path(featuresFile)
  output:
//    set s, file {"*_sorted.bam*"} into bamOutput
//    tuple val(s), file("*_sorted.bam"), file("*_sorted.bam.bai"), file(r) into bamOutput
//    tuple s, file{"${s}*_markup.txt"} into markup
    tuple s, path{"${s}*_annotation.txt"} into annotation
//    tuple s, path{"${s}*_annotation_strings.txt"} into annotationStrings
    tuple s, path{r}, path{"${s}*_features.fasta"} into features mode flatten
    
  script:
    // todo: modularize this chunk
    def bamName = r.name.replaceFirst(".fastq.gz", "").replaceFirst(".fasta.gz", "").replaceFirst(".fastq", "").replaceFirst(".fasta", "")
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
    cut -f2 ${bamName}_annotation.txt | sort | uniq -c > ${bamName}_annotation_strings.txt

    """
} // structure

/*
 * Annotates the collection of a single feature using blat.
 * Output is id_locus.psi file for the input feature sequences.
 */
process blat {
    if(params.nocontainer == "null") { 
        container = "hiroko/blat-hg19:latest"
    }
    //publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  input:
    set s, file(inContig), file(r) from features
    path(refAlleleDir)
  output:
    tuple s, path(inContig), path(r), path{"${s}*.psi"} optional true into psi

  script:
    def nameNoExt = r.name.replaceAll("_features.fasta", "")
    def intervening = nameNoExt.contains("3DP1-2DL4") // skip this one
    def i = nameNoExt.lastIndexOf('_')
    def locus = nameNoExt[i+1..-1]

    """
    if [ "${intervening}" == "false" ]; then
        echo "${intervening}"
        echo "blat: ${r} ${refAlleleDir}/KIR${locus}_nuc.fasta ${s}.psi"
        blat ${r} ${refAlleleDir}/KIR${locus}_nuc.fasta ${s}_${locus}.psi
    fi
    """
} // blat

/*
 * Annotates the blat-annotated sequences into hints.
 * Output is id_locus.gff file for the input feature sequences.
 */
process hints {
    if(params.nocontainer == "null") { 
        container = "chrishah/premaker-plus:18"
    }
    //publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  input:
    set s, file(inContig), file(r), file(p) from psi
  output:
    tuple s, path(inContig), path(r), path{"${s}_*.gff" } optional true into gff

  script:
    def nameNoExt = p.baseName
    def i = nameNoExt.lastIndexOf('_')
    def locus = nameNoExt[i+1..-1]
    """
    echo "blat2hints --in=${p} --out=${s}_${locus}.gff"
    blat2hints.pl --in=${p} --out=${s}_${locus}.gff
    """
} // hints

/*
 * Uses the hints file to annotate with AUGUSTUS gene annotation per locus.
 * Output is *_augustus.gff.
 * requires ENV AUGUSTUS_CONFIG_PATH env to be set
 */
process augustus {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    validExitStatus 0,1
    //publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  input:
    set s, file(inContig), file(r), file(g) from gff
    path(refAlleleDir)    
  output:
    tuple s, path(inContig), path(r), path("*_augustus.gff") optional true into augustus

  script:
    def nameNoExt = g.baseName
    def i = nameNoExt.lastIndexOf('_')
    def locus = nameNoExt[i+1..-1]
    def fullLocus = "KIR" + locus
    """
    #echo "augustus ${s} ${g}"
    #echo $PATH
    /root/augustus/bin/augustus --species=human --UTR=on --strand=both --sample=100 --keep_viterbi=true --alternatives-from-sampling=false --genemodel=partial --hintsfile=${g} --extrinsicCfgFile=extrinsic.ME.cfg --protein=on --introns=on --start=on --stop=on --cds=on --codingseq=on --alternatives-from-evidence=true --proteinprofile=genes/${fullLocus}_prot_msa.prfl ${r} --outfile=${s}_${locus}_augustus.gtf 2> ${s}_augustus_err.txt
    ret=1
    if grep -q gene_id ${s}_${locus}_augustus.gtf; then
        /root/augustus/scripts/gtf2gff.pl < ${s}_${locus}_augustus.gtf --out=${s}_${locus}_augustus.gff --gff3
    else
        exit 0    # not working (todo)
    fi
    """
} // augustus

/*
 * Takes the augustus feature annotation, annotates with IPD-KIR, and
 * makes it suitable for table2asn_GFF.
 * Output is a new GFF and gl.txt
 */
// makes the feature table and annotates wrt IPD-KIR
process alleles {
    if(params.nocontainer == "null") { 
        container = params.container
    }
//    publishDir params.output, mode: 'copy', overwrite: true
    tag { s }

  input:
    set s, file(inContig), file(r), file(g) from augustus // r is fasta, g is gff
    path(refAlleleDir)
  output:
    path("*.gl.txt") optional true into gl
    path("*.ft.txt") optional true into ft
    path("*_mod.gff") optional true into gffMod
  script:
    def nameNoExt = g.name.replaceAll("_augustus.gff", "")
    def i = nameNoExt.lastIndexOf('_')
    def locus = nameNoExt[i+1..-1]
    def fullLocus = "KIR" + locus
    """
    #echo "interpret alleles ${s} ${g}"
    # cdsexons
    /root/augustus/scripts/getAnnoFasta.pl --seqfile ${r} ${g} 2> ${s}_${locus}_getAnnoFasta_err.txt
    # aa and codingseq
    /root/augustus/scripts/getAnnoFastaFromJoingenes.py -g ${r} -3 ${g} -o ${s}_${locus}_augustus 2> ${s}_${locus}_getAnnoFastaFromJoingenes_err.txt
    gff2ftGene.groovy -g ${fullLocus} -i ${refAlleleDir} -f ${g} -s ${r} -o ${s}_${fullLocus}.ft.txt -l ${s}_${fullLocus}.gl.txt 2> ${s}_${locus}.alleles_err.txt
    augustus2NCBI.groovy -g ${fullLocus} -i ${refAlleleDir} -f ${s}_${locus}_augustus.gff -s ${r} -l ${s}_${fullLocus}.gl.txt 2> ${s}_${locus}_aug2NCBI_err.txt
    """
} // alleles

allModGFFs = gffMod.collectFile(name: 'mod.gff', newLine: true)
modGFFwFasta = allModGFFs.combine(annotateReads)

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

/*
 * Combine the per-feature GFFs into the per contig annotation.
 * Run through table2asn_GFF to create the NCBI gff, sqn, and a report.
 */
process combineGFF {
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir params.output, mode: 'copy', overwrite: true
    input:
        tuple path(modGFF), val(desc), path(inContig) from modGFFwFasta
        file sbtF from params.sbt
    output:
        path("*.sqn")
        path("*.gbf")
        path("*.dr")
        path("*.tbl")
        path("*.err")
    script:
    """
    echo "${modGFF}, ${desc}, ${inContig}"
    combineGFF.groovy -i mod.gff -o .
    linux64.table2asn_GFF -augustus-fix -f ${desc}.gff -i ${inContig} -outdir . -genbank -verbose -euk -V b -Z  -t ${sbtF} -j "[organism=Homo sapiens]"
    gbf2tbl.pl ${desc}.gbf
    """
} // combineGFF

// combineGLs todo

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf('.'))

}
