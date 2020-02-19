#!/usr/bin/env groovy

/*
 * alignment2ProbePairs
 * 
 * For each read in each sam file in the input path, 
 * record the order & location of the probe pairs.
 * Represent each contig as an ordered list of probe pairs; the probe names
 * are separated by '-'. Each probe pair
 * is represented by a single character from a 32 character set.
 * Then, can compare via longest common substrings, etc.
 * The input is a directory of *sorted* SAM files.
 * 
 * e.g., alignment2ProbePairs.groovy -d . -o 15_markerHap_v1.txt
 * 
 * Requires
 *  guava.jar: https://github.com/google/guava
 *    http://google.github.io/guava/releases/19.0/api/docs/com/google/common/collect/Table.html
 *
 * @author Dave Roe
 * @todo pass in the fasta regex?
 * @todo normalize wrt direction; MHC, not CHM
 *
 * bowtie2 -a --end-to-end
 * samtools view -b -S contgs1.sam > contgs1.bam
 * samtools sort contgs1.bam -O SAM -o contgs1_sorted.sam
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// full 52 character set (from ASCII?)
//accessionList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/']
// 32+ from the 52 character set (from ASCII?)
accessionList = ['A', 'N', 'O', 'P', 'Q', 'e', 'f', 'g', 'h', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/']
// use these to normalize the direction of the markup
// same direction as the sequences in GeneBank, etc.
// todo: don't hard code this (how?)
//directionCluesList = ["MHC", "FGH", "SCT", "WXY", "JKL"]
// todo: make all combinations?
directionCluesList = ["MH", "HC", "FG", "GH", "SC", "CT", "WX", "XY", "JK", "KL", "CD",
                      "LF", "LM", "Lb", "DE", "Dc", "EL", "cL", "UV", "VW", "XY", "CA",
                      "AN", "JO", "IJ", "JK", "KL", "CI", "RL", "IR", "PQ"]

// map probe Name -> accession number stored as Byte
// probe pair name is <probe>-<probe> (todo: change this to a Table)
probeAccMap = ['13-4':'H', '3-9':'X', '12-2':'R', '7-13':'M', '10-2':'K', '2-4':'k', '3-11':'D', '7-12':'b', '6-7':'a', '2-15':'i', '11-7':'E', '1-1':'Z', '6-4':'B', '2-7':'L', '3-12':'I', '1-13':'G', '3-14':'T', '12-10':'J', '14-5':'U', '7-15':'m', '5-8':'V', '11-2':'c', '5-4':'l', '7-1':'F', '4-3':'C', '5-6':'d', '8-3':'W', '5-3':'j', '9-10':'Y', '7-4':'S']
// map probe Name -> accession number stored as Byte
// use this when not predefined
//HashMap<String, Character> probeAccMap = new HashMap()

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3

// thing that probably won't change per run
maxProbeDistance = 10000   // probe pairs must be this distance or less
minProbePairSize = 1 // only output strings with this many characters
err = System.err

OptionAccessor options = handleArgs(args)
boolean outputMatrix = true
err.println "options.m=" + options.m //todo
if(options.m == "0") {//todo: this isn't working
	outputMatrix = false
}

// open output file
outWriter = new PrintWriter(new File(options.o).newOutputStream(), true)

// loop through each sam file within the input directory
dir = new File(options.d)
fp = ~/.*.sam/
err.println "looking for " + fp
// loop through the sam files
dir.eachFileMatch(fp) { inFile ->
	err.println "processing ${inFile}"
	// convert sam file to marker-pair string
	processSam(inFile, probeAccMap, directionCluesList, outWriter)
} // each input file
outWriter.println ""
err.println probeAccMap.keySet().size() + " probe pairs"

// output the probe accession mapping, for fyi and
// so they can be reused
if(outputMatrix == true) {
	outWriter.println probeAccMap.keySet().size() + " probe pairs"
	outWriter.print "["
	first = true
	probeAccMap.each { probePair, probeAcc ->
		if(first) {
			first = false
		} else { 
			outWriter.print ", "
		}
		outWriter.print "'${probePair}':'${probeAcc}'"
	} // each probe and accession
	outWriter.println "]"
}

/*
 * processSam
 * Prints the probe pairs sequentially, representing them with
 * the accessions in the map. It populates the probe to accession map.
 * 
 */
void processSam(File samFile, HashMap<String,Character> probeAccMap,
				ArrayList directionCluesList, PrintWriter outWriter) {
	reader = new FileReader(samFile)
	// row, column, cell: read, location, probe name
	HashBasedTable<String, Integer, String> probeHitTable = HashBasedTable.create()
	// populate probeHitTable
	reader.eachLine { line ->
        ArrayList cols = line.split('\t')
        String probe = cols[0].trim() // probe name
		 // skip blank or header rows
		String read = cols[2].trim() // read/contig name
        if((probe == null) || (probe == "") || (probe.startsWith("@")) ||
		   (read == "*")) {
            return
        }
		Integer location = cols[3].toInteger()

		probeHitTable.put(read, location, probe)
	} // each line of sam

    annotateReads(samFile, probeHitTable, probeAccMap, directionCluesList,
				  outWriter)
} // processSam

/*
 * annotateReads
 * @param probeHitTable row, column, cell: read, location, probe name
 */ 
void annotateReads(File samFile,
				   HashBasedTable<String, Integer, String> probeHitTable,
				   HashMap<String,Character> probeAccMap,
				   ArrayList directionCluesList, 
				   PrintWriter outWriter) {

    Long distanceSum = 0
    Integer distanceCount = 0
    Integer distanceMax = 0

	// annotate each read and populate probeAccMap
	probeHitTable.rowKeySet().each { read ->
		if(debugging <= 2){
			err.println "read ${read}"
		}
		Map<Integer, String> colMap = probeHitTable.row(read)
		locationSet = colMap.keySet().sort() as Set
		Integer previousLocation = 0
		String previousProbe = ""
		Integer location = null
		String outString1 = ""
		locationSet.each { inLoc ->
			location = new Integer(inLoc)
			if(debugging <= 2) { 
				err.println "location ${location}"
				err.println "previousLocation ${previousLocation}"
			}
			Integer distance = location - previousLocation
            distanceSum += distance
            distanceCount++
            //if((distance > distanceMax) && (previousLocation > 30000)) {
            if((distance > distanceMax) && (previousLocation != 0)) {
                err.println "found max distance ${distance} (${previousLocation}-${location}) in ${read}"
                distanceMax = distance
            }
			boolean unexpectedGap = (previousLocation == 0) ?
				false : (distance > maxProbeDistance)
			if((previousLocation == 0) || unexpectedGap) {
				if(unexpectedGap) {
					err.println "WARNING: unexpected gap (${distance} bp) in ${read} (${samFile}) between ${previousLocation} (${colMap[previousLocation]}) and ${location} (${colMap[location]})"
				}
				previousLocation = location
				previousProbe = colMap[location]
				return
			}
			probe = colMap[location]
			
			probePair = previousProbe + "-" + probe
			Character probeAcc = getProbePairAcc(probeAccMap, probePair)
			if(probeAcc == null) {
				probeAcc = accessionList.pop()
				probeAccMap[probePair] = probeAcc
				if(debugging <= 2) { 
					err.println "$location $probePair = $probeAcc"
					err.println "adding to accessionList: " + accessionList.join(",")
				}
			}
			if(debugging <= 2) {
				err.println "${probePair}=${probeAcc}"
			}
			outString1 += probeAcc

			previousProbe = probe
			previousLocation = location
		} // each location

		if((outString1.length()+1) != locationSet.size()) {
			err.println "ERROR: markup size doesn't match locations for ${read}"
			err.println "annotateReads: in outString1 (length ${outString1.length()}) = ${outString1}"
			err.println "extractDNAFeaturesFromAll: locations (length ${locationSet.size()}) =" + locationSet
			//System.exit(1)
		}
		outputAnnotation(read, outString1, locationSet,
						 directionCluesList, outWriter)
	} // each read
    avg = distanceCount ? distanceSum/distanceCount : 0
    err.println "${distanceCount} reads; average=${avg}; max=${distanceMax}"
} // annoateReads


/* 
 * outputAnnotation
 * 
 */
void outputAnnotation(String read, String outString1,
					  Set locationSet, ArrayList directionCluesList,
					  PrintWriter outWriter) {
	if(outString1.length() >= minProbePairSize) {
		// standardize the direction of the markup
		reverseIt = true
		directionCluesList.each { clue ->
			if(outString1.contains(clue)) {
				reverseIt = false
				return
			}
		} // each clue
		if(reverseIt) { 
			outString1 = outString1.reverse()
		}
		outWriter.print "${read}\t${outString1}\t"
		outWriter.println locationSet.sort().join(',')
	}
} // outputAnnotation


/* 
 * Returns an accession character given a probe pair
 * separated by a '-'.
 */
Character getProbePairAcc(HashMap<String,Character> probeAccMap, String pp) {
	Character ret = probeAccMap[pp]
	if(debugging <= 2){
		err.println "getProbePairAcc: ${pp} = ${ret}"
	}

	if(ret == null) {
		(p1, p2) = pp.split("-")
		ret = probeAccMap[new String(p2 + "-" + p1)]
	}
	return ret
} // getProbePairAcc

/*
 * handleArgs

 * 
 * Parse and return the input.
 *
 * @todo e.g., here
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'alignment2ProbePairs.groovy [options] ', header:'Options:')
    cli.d(longOpt:'directory', args:1, argName:'dir', 'directory with SAM files',
		  required: true)
    cli.m(longOpt:'matrix', args:1, argName:'matrix', 'output the matrix',
		  required: false)
    cli.o(longOpt:'output', args:1, argName:'out',
		  'output file containing marker strings', required: true)
    OptionAccessor options = cli.parse(args)

    return options
} // handleArgs
