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
params.home = baseDir
home = params.home + "/"
params.raw = home + "/raw/"
raw = params.raw + "/"
params.output = home + "/output"
output = params.output + "/"
params.off = 0
fqNameSuffix = "fastq.gz"          // extension on the file name (todo: expand this)
params.container = "droeatumn/kass:latest"
params.nocontainer = "null"
params.canuPB = "-pacbio-hifi"
params.threads = "8"

// things that probably won"t change per run
fqPath = raw + "/*" + fqNameSuffix
capFile = file("${home}/input/cap.fasta") // capture probes
markerCapFile = file("${home}/input/markers_wCap.fasta") // gene markers + capture probes
featuresFile = file("${home}/input/features.txt") // markup features
haps = home + "${home}/input/HapSet23_v1.txt"
alignProbesFile = file("${home}/src/alignment2ProbePairs.groovy")
annotateFile = file("${home}/src/annotateMarkup.groovy")

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
    if(params.nocontainer == "null") { 
        container = params.container
    }
    publishDir output, mode: 'copy', overwrite: true
//doesn't work    publishDir "*_off-kir.fastq.gz", mode: 'copy', overwrite: true
//doesn't work    publishDir "*_kir.fastq", mode: 'copy', overwrite: true
    input:
        tuple s, path(fa) from fqs
        path(markerCapFile)
    output:
	    tuple s, file{"*_kir.fastq"} into kirFastqs
        tuple s, file{"*_off-kir.fastq.gz"} into offkirFastqs optional true

    script:
        def offFile = ""
        def offStr = ""
        if(params.off != 0) {
            offFile = s + "_off-kir.fastq"
            offStr = "out="  + offFile
        }
    """
    bbduk.sh in=${fa} ${offStr} outm=${s}_kir.fastq ref=${markerCapFile} k=25 maskmiddle=f overwrite=t rename=t nzo=t rcomp=t ignorebadquality=t
    # remove empty files
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
  if(params.nocontainer == "null") { 
      container = params.container
  }
  //publishDir output, mode: 'copy', overwrite: true
  input:
    tuple s, path(fq) from kirFastqs
  output:
	tuple s, path{"${s}*.fasta"} into correctedReads
//    	tuple s, path{"*.fasta.gz"} into correctedReads mode flatten
	
    """
    lorma.sh ${fq}
    mv final.fasta ${s}-corrected.fasta

    binBBFasta.groovy -i ${s}-corrected.fasta -o .

#    gzip ${s}-corrected*.fasta
    """
} // correct

/*
 * assemble
 * 
 * @todo change genome size for genes
 * @todo rename the contig names to eliminate duplications
 */
process assemble {
  if(params.nocontainer == "null") { 
    container = params.container
  }
  publishDir output, mode: 'copy', overwrite: true //todo
  errorStrategy 'ignore'
    
  input:
    tuple s, path(cr) from correctedReads
  output:
    tuple s, path{"*.contigs.fasta.gz"} into assembly

    script:
    //def s2 = cr.name.replaceFirst("-corrected", "").replaceFirst(".fasta.gz", "")
    """
    id=""
    firstID=""
    FILES="*-corrected*.fasta"
    for bFile in \$FILES; do
        echo \$bFile
        id=\$(basename \$bFile)
        # '%' Means start to remove after the next character;
        id=\${id/-corrected/}
        id=\${id%.fasta}
        echo \$id
        if [ "\$firstID" == "" ]; then
            firstID=\$id
        fi
        canu -p \$id -d \$id genomeSize=200k ${params.canuPB} \$bFile || true
        cp \$id/\$id.contigs.fasta . || true
        deep.pl replace '>tig' ">\${id}_tig" "\${id}.contigs.fasta"
    done
    cat *.contigs.fasta > tmp.fasta
    orient.groovy -i tmp.fasta -p ${capFile} -o tmp2.fasta
    reformat.sh in=tmp2.fasta out=\$firstID.contigs.fasta fastawrap=1000000 overwrite=true
    gzip \$firstID.contigs.fasta
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

