#!/usr/bin/env groovy
/*
 * Make some modifications to the table (tbl) file.
 *
 * @author Dave Roe
 * @todo add "<" and ">" for partial; requires a refactor to process by section
 */
import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
import groovy.transform.Field
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3

// thing that probably won't change per run
err = System.err

@Field final OptionAccessor options = handleArgs(args)

// description (as in the fasta), position -> gff feature
HashBasedTable<String, Integer, String> idPosTblTable = HashBasedTable.create()
inFileName = options.i
FileReader iReader = new FileReader(inFileName)
PrintWriter writer = new PrintWriter(new File(options.o).newOutputStream(), true)

loadTbl(iReader, writer)

void loadTbl(tReader, PrintWriter writer) {
    l1 = tReader.readLine()
    l1Split = l1.split('\t')
    start = 0 // of gene
    end = 0
    partial = false
    currentSection = null  // 0: gene, 1: mRNA, 2: CDS
    gene = null
    for (String l2 = tReader.readLine(); l2 != null; l2 = tReader.readLine()) {
        l2Split = l2.split('\t')
        l2SplitSize = l2Split.size()
        att = ""
        if(l2SplitSize > 3) {
            att = l2Split[3..-1].join('\t')
        }
        if(debugging <= 2) {
            err.println "loadTbl: l2=${l2}"
            err.println "loadTbl: l2Split size=${l2Split.size()}"
            err.println "loadTbl: att=${att}"
        }
        if(l1Split[0].startsWith(">Feature")) { // e.g., >Feature |MN167504_3DL2_160218-171148|
            (f, name) = l1Split[0].split(' ')
            l1Split[0] = f + " |" + name + "|"
        }
        if((l2SplitSize > 3) && att.contains("partial")) {
            partial = true
        } else if((l2SplitSize > 3) && att.startsWith("gene") &&
                  (currentSection == 0)) {
            (geneStr, gene) = att.split('\t')
        } else if((l2SplitSize > 3) && att.startsWith("gene") &&
                  (currentSection == 1)) {
            // add 3' UTR unless partial 3DL2; previous line
            if(!((partial == true) && gene.contains("3DL2"))) {
                l1Split[1] = end // add 3' UTR
            } else if(gene.contains("3DL2")) {
                start = l1Split[0].toInteger()
                end = l1Split[1].toInteger()
                err.println "in 1 3DL2 ${end} - ${start} - 1: ${end-start+1}"//todo
                err.println "l1Split=${l1Split}"//todo
                if(end-start+1 == 329) {
                    err.println "changing 3DL2 ${end-start+1} into 294"//todo
                    l1Split[1] = start + 294 - 1
                }
            }
            continue; // skip allele in mRNA
        } else if((l2SplitSize > 3) &&
                  (att.startsWith("allele") || att.startsWith("gene")) &&
                  ((currentSection == 1) || (currentSection == 2))) {
            if(gene.contains("3DL2")) {
                start = l1Split[0].toInteger()
                end = l1Split[1].toInteger()
                err.println "in 1 or two 3DL2 ${end} - ${start} - 1: ${end-start+1}"//todo
                err.println "l1Split=${l1Split}"//todo
                if((end-start+1 == 329) || (end-start+1 == 625)) {
                    err.println "changing 3DL2 ${end-start+1} into 294"//todo
                    l1Split[1] = start + 294 - 1
                }
            }
            continue; // skip gene and allele in CDS
        } /*else if((l2SplitSize > 3) && att.startsWith("codon_start") &&
                  (currentSection == 2)) {
            (cs, num) = att.split('\t')
            if(num != "1") { 
                // output the previous line 
                writer.println l1Split.join('\t')
                // make l1 the codon_start line
                l1Split = l2Split.collect()
            }
        } */ else if((l2SplitSize > 3) && att.startsWith("product") &&
                  (currentSection == 2)) {            
            if(gene.contains("DP")) { // skip CDS in pseudo genes 2DP1, 3DP1
                continue;
            }
            cs = l1Split[3]
            num = l1Split[4]
            if(cs.startsWith("codon") && (num != "1")) { 
                writer.println l2Split.join('\t') // print the product line first
            } else { 
                l1Split = l2Split.collect() // product is now l1
            }
            // add transl_table
            l2Split[3] = "transl_table"
            l2Split[4] = "1"
        } else if((l2SplitSize > 3) && (gene.contains("DP") == true) &&
                  (currentSection == 2)) { // skip CDS stuff in pseudos 2DP1, 3DP1
            continue;
        } else if((l2SplitSize > 2) && (l2Split[2] == "gene")) {
            start = l2Split[0].toInteger()
            end = l2Split[1].toInteger()
            partial = false
            currentSection = 0
            gene = null
        } else if((l2SplitSize > 2) && (l2Split[2] == "mRNA")) {
            // add 5' UTR unless partial 3DL3
            if(!((partial == true) && gene.contains("3DL3"))) {
                l2Split[0] = start
            }
            currentSection = 1
        } else if((l2SplitSize > 2) && (l2Split[2] == "CDS")) {
            currentSection = 2
            if(gene.contains("DP")) { // skip CDS in pseudo genes 2DP1, 3DP1
                continue; 
            }
        }

        writer.println l1Split.join('\t')
        l1Split = l2Split
    } // each line
    writer.println l1Split.join('\t')

} // loadTbl

/*
 * handleArgs

 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'modifyTbl.groovy [options] ', header:'Options:')
    cli.i(longOpt:'input', args:1, argName:'in', 'tbl input file', required: true)
    cli.o(longOpt:'output', args:1, argName:'out', 'tbl output file', required: true)
    OptionAccessor options = cli.parse(args)

    return options
} // handleArgs
