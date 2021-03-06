#!/usr/bin/env groovy
/*
 * Combine multi-contig, multi-gene feature table files into one per contig.
 *
 * Removes db_xref lines.
 *
 * @author Dave Roe
 * @todo skip non-feature stuff like source, misc_feature, etc
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

// overall goal is to consolidate feature table per id (description)
// description (as in the fasta), position -> feature table
HashBasedTable<String, Integer, String> idPosFtTable = HashBasedTable.create()
FileReader iReader = new FileReader(options.i)

loadFT(iReader, idPosFtTable)

ArrayList<String> ids = idPosFtTable.rowKeySet()
ArrayList<Integer> positions = idPosFtTable.columnKeySet().sort()
if(debugging <= 3) {
	err.println "${ids.size()} ids: ${ids.join(',')}"
    err.println "${positions.size()} positions: ${positions.join(',')}"
}
ids.each { i ->
    String featuresFull = new String()
    positions.each { pos ->
        featureTmp = idPosFtTable.get(i, pos)
        if((featureTmp != null) && (featureTmp.trim() != "")) {
            featuresFull += featureTmp
        }
    }
    outName = i
    PrintWriter writer = new PrintWriter(new File("${outName}.ft.txt").newOutputStream(), true)
    writer.println ">Feature ${i}"  // e.g., >Feature MN167504
    writer.println featuresFull
    writer.close()
} // each contig

// load feature table into the hash table
void loadFT(tReader, idPosFtTable) {
    String desc = null
    Integer position = null
    String featureString = ""
    String geneLine = ""
    Boolean pastGene = false
    //    while(l = tReader.readLine()) {
    for (String l = tReader.readLine(); l != null; l = tReader.readLine()) {
        lSplit = l.split('\t')
        if(debugging <= 2) {
            err.println "loadFT: l=${l}"
            err.println "loadFT: lSplit size=${lSplit.size()}"
        }
        
        if(l.startsWith(">Feature")) {
            if(desc != null) {
                if(debugging <= 3) {
	                err.println "(${desc}, ${position}) = ${featureString}"
                    err.println "geneLine = ${geneLine}"
                }
                featureString = geneLine + '\n' + featureString
                idPosFtTable.put(desc, position, featureString)
                desc = null
                position = null
                featureString = ""
                geneLine = ""
            }
            (fStr,id) = l.split('Feature ')
            // separation location
            idx = id.lastIndexOf('_')
            desc = id[0..idx-1]
            // remove locus from description
            lIdx = desc.lastIndexOf('_')
            desc = desc[0..lIdx-1]
            
            posTmp = id[idx+1..-1]
            (posTmp, end) = posTmp.split('-')
            position = posTmp.toInteger()
            if(debugging <= 3) {
                err.println "loadFT: setting postion to ${position}"
            }
            pastGene = false
            continue
        } else if(l.contains("db_xref")) { // remove db_xref
            continue
        } else if(desc == null) {  // skip stuff before first feature; source, etc.
            // todo: this only works for the first one (todo)
            continue
        } else if(l.contains("source")) {    // skip stuff before first feature; source, etc.
            // this is temporary hard-coding (todo) *
            // skip the next 4 lines
            (1..7).each {
                tReader.readLine()
            }
            continue
        } else if((lSplit.size() > 2) && (lSplit[2] == "gene")) {
            if(debugging <= 3) {
                err.println "gene line found: ${l}"
            }
            geneLine = l
            continue
        } else if(l.contains("mRNA")) {
            pastGene = true
        } else if(l.contains("\tallele\t") && (pastGene == true)) {
            // remove the 'allele' from mRNA and CDS sections
            return
        }

        featureString += l + '\n'
    } // while
    if(debugging <= 3) {
	    err.println "(${desc}, ${position}) = ${featureString}"
        err.println "geneLine = ${geneLine}"
    }
    featureString = geneLine + '\n' + featureString
    idPosFtTable.put(desc, position, featureString)
} // loadFT

/*
 * handleArgs

 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'combineFT.groovy [options] ', header:'Options:')
    cli.i(longOpt:'input', args:1, argName:'in', 'feature table input file', required: true)
    //    cli.f(longOpt:'fasta', args:1, argName:'fa', 'fasta input file', required: true)
//    cli.o(longOpt:'output', args:1, argName:'out', 'output file', required: true)
    OptionAccessor options = cli.parse(args)

    return options
} // handleArgs
