#!/usr/bin/env groovy
/*
 * Combine multi-contig, multi-gene GFF into one per contig.
 *
 * @author Dave Roe
 * @todo could take a fasta file and only output the gff for that one
 */
import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
import groovy.transform.Field
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3

// thing that probably won't change per run
err = System.err

@Field final OptionAccessor options = handleArgs(args)

// overall goal is to consolidate gffs per id (description)
// description (as in the fasta), position -> gff feature
HashBasedTable<String, Integer, String> idPosGffTable = HashBasedTable.create()
FileReader iReader = new FileReader(options.i)

loadGFF(iReader, idPosGffTable)

ArrayList<String> ids = idPosGffTable.rowKeySet()
ArrayList<Integer> positions = idPosGffTable.columnKeySet().sort()
if(debugging <= 3) {
	err.println "${ids.size()} ids, ${positions.size} positions"
}
ids.each { i ->
    String featuresFull = new String()
    positions.each { pos ->
        featureTmp = idPosGffTable.get(i, pos)
        if((featureTmp != null) && (featureTmp.trim() != "")) {
            featuresFull += featureTmp
        }
    }
    PrintWriter writer = new PrintWriter(new File("${options.o}${File.separator}${i}.gff").newOutputStream(), true)
    writer.println featuresFull
    writer.close()
} // each contig

// load GFF into the table
void loadGFF(iReader, idPosGffTable) {
    String desc = null
    Integer position = null
    String featureString = ""
    iReader.eachLine { l ->
        if(l.startsWith("# Feature")) {
            if(desc != null) {
                if(debugging <= 3) {
	                err.println "(${desc}, ${position}) = ${featureString}"
                }
                idPosGffTable.put(desc, position, featureString)
                desc = null
                position = null
                featureString = ""
            }
            (fStr,id) = l.split('Feature:')
            (desc, posTmp) = id.split('-')
            position = posTmp.toInteger()
            return
        } else if(l.contains("\tlocus_tag\t")) { // remove locus_tag
            return
        } /*else if(l.contains("\tprotein_id\t")) { // no protein_id in tbl (ncbi)
            return*/
        featureString += l + '\n'
    } // each input line
    if(debugging <= 3) {
	    err.println "(${desc}, ${position}) = ${featureString}"
    }
    idPosGffTable.put(desc, position, featureString)
} // loadGFF

/*
 * handleArgs

 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'combineGFF.groovy [options] ', header:'Options:')
    cli.i(longOpt:'input', args:1, argName:'in', 'GFF input file', required: true)
    //    cli.f(longOpt:'fasta', args:1, argName:'fa', 'fasta input file', required: true)
    cli.o(longOpt:'output', args:1, argName:'out', 'output directory', required: true)
    OptionAccessor options = cli.parse(args)

    return options
} // handleArgs
