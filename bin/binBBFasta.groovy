#!//usr/local/sdkman/candidates/groovy/current/bin/groovy
//#!/usr/bin/env groovy

/*
 * binBBFasta
 *
 * Given a FASTA file with bbduk-labeled header lines, bin the reads into 
 * separate files based on the KPI gene markers in the header.
 *
 * binBBFasta.groovy -i <input fasta uncompressed>  <output location for fasta uncompressed files>
 * e.g, binBBFasta.groovy -i Kir5_P1_8_999F_corrected_kir.fasta -o .
 *
 * @author Dave Roe
 */

import groovy.io.*
import groovy.cli.commons.OptionAccessor
import groovy.cli.commons.CliBuilder

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3

// things that probably won't change per run
err = System.err
separatorCharacter = '-'  // separates the probe pairs
newline = System.getProperty('line.separator')

OptionAccessor options = handleArgs(args)
err.println options//todo
binReads(options.i, options.o)
// end main

/*
 * binReads
 *
 * Bin the input fasta lines into their KIR gene regions.
 * 
 * @param faFileName input fasta file
 * @param outDirName location to place the output
 */
def void binReads(String faFileName, String outDirName) {
    if(debugging <= 3) {
        err.println "binReads(faFileName=${faFileName}, outDirName=${outDirName})"
    }
	counter = 0
	progressLevel = 1000 // show progress every x markers

    // open input file with fasta
    FileReader fastaReader = new FileReader(new File(faFileName))
	String outID = outDirName + "/" +
		faFileName.split("/")[-1].replaceFirst(".fasta", "").replaceFirst(".fa", "")
	// map: regions -> PrintWriter for per-region output files
	HashMap<String,PrintWriter> regionWriterMap = new HashMap()

	// assigned region from each read's header
	HashSet<String> headerRegions = new HashSet()
	Queue<String> currentTwoLines = new LinkedList() as Queue
    fastaReader.eachLine { line ->
		if(line.startsWith(">")) { //todo: parameterize this
			writeTwoLines(outID, regionWriterMap, headerRegions, currentTwoLines)
			headerRegions.clear()
			currentTwoLines.clear()
			
			(header, headerRegions) = parseHeaderCombinedCount(line)
			if(debugging <= 2) { 
				err.println "${headerRegions.size()} regions for $header"
			}
			currentTwoLines.add(header)
		} else {
			currentTwoLines.add(line)
		} // header vs other three lines
    	if((++counter % progressLevel) == 0) {  // print progress '.'
	    	err.print "."
    	}
	} // each fasta line
	writeTwoLines(outID, regionWriterMap, headerRegions, currentTwoLines)
	regionWriterMap.values().each { w ->
		w.close()
	}
	fastaReader.close()
	
    if(debugging <= 1) {
        err.println "binReads: return"
    }
    return
} // binReads

/*
 * Write the two fasta lines to each region.
 * The file name convention is outDir/${id}_${region}.fasta
 * 
 * @param id String containing the path plus the id of the file name.
 */
void writeTwoLines(String id, HashMap<String, PrintWriter> regionWriterMap,
					HashSet<String> headerRegions,
					Queue<String> currentTwoLines) { 
	PrintWriter outWriter
	headerRegions.each { region ->
		outWriter = regionWriterMap[region]
		if(outWriter == null) {
			name = "${id}_${region}.fasta"
			outWriter = new PrintWriter(new File(name).newOutputStream(), true)
			regionWriterMap[region] = outWriter
		}
		currentTwoLines.each { line ->
			outWriter.println line
		}
	}
} // writeTwoLines

/* 
 * parseHeaderCombinedCount
 *
 * Return the normal header plus a HashSet of the regions to which to 
 * assign this read.
 * 
 * e.g., @m180626_142523_42289_c101463372550000001823304808281891_s1_p0/49812/ccs	2DL4=1
 * @return List with a String of the original header and a HashSet<String> of the regions with the most hits
 */
ArrayList parseHeaderCombinedCount(String line) {
	if(debugging <= 1) {
		err.println "parseHeaderCombinedCount(line=$line)"
	}
	TreeSet<String> regionSet = new TreeSet()
	TreeSet<String> hitNameSet = new TreeSet()
	HashMap<String, Integer> headerCounts = new HashMap()

	ts = line.split('\t')
	retHeader = ts[0] //the original header
	maxCount = 0
	if(ts.size() > 1) { // if there was at least one label
		ts[1..-1].each { t ->
			if(!t.contains("=")) {
				return
			}
			(region, rest) = t.split("=")
            if(region =~ /^[23]D/) { // check for genes
			    currentVal = headerCounts[region] ?: 0
			    headerCounts[region] = ++currentVal
			    if(currentVal > maxCount) {
				    maxCount = currentVal
			    }

			    regionSet.add(region)
            }
		}
	}

	if(debugging <= 1) {
		err.println "parseHeaderCombinedCount: return " + regionSet
	}
	return [retHeader, regionSet]
} // parseHeaderCombinedCount

/* 
 * parseHeaderMaxCount
 *
 * Return the normal header plus a HashSet of the regions to which to 
 * assign this read.
 * 
 * e.g., @m180626_142523_42289_c101463372550000001823304808281891_s1_p0/49812/ccs	2DL4=1	2DL4=1	2DL4=1	2DL4=1	2DL4=1	2DS1=1	2DP1-2DS1=1	2DL2-3DP1=1	2DL4=1	2DL3-2DL5B=1	2DS5-2DS1=1	2DS5-2DP1=1	2DS1=1	2DL4-3DS1=1	2DL3-2DP1=1	2DS2-2DL2=1
 * @return List with a String of the original header and a HashSet<String> of the regions with the most hits
 */
ArrayList parseHeaderMaxCount(String line) {
	if(debugging <= 1) {
		err.println "parseHeaderMaxCount(line=$line)"
	}
	HashSet<String> headerRegions = new HashSet()
	HashMap<String, Integer> headerCounts = new HashMap()
	
	ts = line.split('\t')
	retHeader = ts[0] //the original header
	maxCount = 0
	if(ts.size() > 1) { // if there was at least one label
		ts[1..-1].each { t ->
			if(!t.contains("=")) {
				return
			}
			(region, rest) = t.split("=")
			currentVal = headerCounts[region] ?: 0
			headerCounts[region] = ++currentVal
			if(currentVal > maxCount) {
				maxCount = currentVal
			}
		}
	}
	headerCounts.each { region, count ->
		if(count == maxCount) {
			if(debugging <= 3) {
				err.println "parseHeaderMaxCount: adding $region for $retHeader"
				err.println "parseHeaderMaxCount: maxCount=$maxCount"
			}

			headerRegions.add(region)
		}
	}

	if(debugging <= 1) {
		err.println "parseHeaderMaxCount: return"
	}
	return [retHeader, headerRegions]
} // parseHeaderMaxCount

/* 
 * parseHeaderAll
 *
 * Return the normal header plus a HashSet of the regions to which to 
 * assign this read.
 * 
 * e.g., @m180626_142523_42289_c101463372550000001823304808281891_s1_p0/49812/ccs	2DL4=1	2DL4=1	2DL4=1	2DL4=1	2DL4=1	2DS1=1	2DP1-2DS1=1	2DL2-3DP1=1	2DL4=1	2DL3-2DL5B=1	2DS5-2DS1=1	2DS5-2DP1=1	2DS1=1	2DL4-3DS1=1	2DL3-2DP1=1	2DS2-2DL2=1
 */
ArrayList parseHeaderAll(String line) {
	HashSet<String> headerRegions = new HashSet()
	if(debugging <= 1) {
		err.println "parseHeaderAll(line=$line)"
	}
	ts = line.split('\t')
	retHeader = ts[0] //the original header
	if(ts.size() > 1) { // if there was at least one label
		ts[1..-1].each { t ->
			if(!t.contains("=")) {
				return
			}
			(region, rest) = t.split("=")
			headerRegions.add(region)
		}
	}

	if(debugging <= 1) {
		err.println "parseHeaderAll: return"
	}
	return [retHeader, headerRegions]
} // parseHeaderAll

/*
 * loadFeatureNomenclatureFile
 * Loads the files containing the markup sequences and their features.
 * 
 * See the top-level documentation for input file.
 * Removes "3DL3~" and "~3DL2" from the markup strings.
 *
 * @return four Maps
 *         markup string -> nomenclature string
 *         integer size -> string markup
 *         region name -> gene names
 *         gene name -> markup
 */
ArrayList<HashMap> loadFeatureNomenclatureFile(String refFileName) {
	if(debugging <= 1) {
		err.println "loadFeatureNomenclatureFile(${refFileName})"
	}

	FileReader reader = new FileReader(refFileName)
	// markup string -> nomenclature string
    HashMap<String,String> markNomenMap = new HashMap()
	// integer size -> list of markup strings
    HashMap<Integer,ArrayList<String>> sizeMarkMap = new HashMap()
	// map of region to gene list; e.g., cA01 -> 3DL3~2DL3~2DP1~2DL1~3DP1
	HashMap<String, String> genesRegionMap = new HashMap()
	// map of (gene) nomenclature to markup
	HashMap<String,TreeSet<String>> nomenMarkupMap = new HashMap()
	
	header = reader.readLine() // header
	if(debugging <= 1) {
		err.println "loadFeatureNomenclatureFile: header=${header}"
	}
	reader.eachLine { line ->
		if(debugging <= 1) {
			err.println "loadFeatureNomenclatureFile: line=${line}"
		}
		if((line == null) || (line == "")) {
			return
		}
        ArrayList cols = line.split('\t')
		err.println "cols=" + cols//todo
        String mark = cols[0].trim() // markup
		//mark = mark.replaceFirst("3DL3~", "").replaceFirst("~3DL2", "")
        String resolution = cols[1].trim() // resolution
		String nomenRegion = cols[2].trim()
        String nomenListStr = cols[3] ? cols[3].trim() : null // gene nomanclature (classic)
		// only use the gene rows, not the regions, etc
		if((resolution == null) || !resolution.contains("gene")) {
			return
		}
		err.println "nomenListStr=${nomenListStr}"//todo
		if(nomenListStr == null) {
			return
		}

		nomenListStr.split("/").each { nomen ->
			// skip the fusions and intergenes
			// e.g., 2DS2\2DS3 and 3DP1-2DL4
			if(nomen.contains("\\") || nomen.contains("-")) {
				return
			}
			markNomenMap[mark] = nomen
			l = mark.length()
			ArrayList<String> mList = sizeMarkMap[l]
			if(mList == null) {
				mList = new ArrayList()
			}
			mList.add(mark)
			sizeMarkMap[l] = mList

			// some genes have markup alleles, not one fixed markup
			nomenMarkupList = nomenMarkupMap[nomen] ? nomenMarkupMap[nomen] : new ArrayList()
			nomenMarkupList.add(mark)
			if(debugging <= 1) { 
				err.println "loadFeatureNomenclatureFile: adding ${nomenMarkupList} for ${nomen} in nomenMarkupList"
			}
			nomenMarkupMap[nomen] = nomenMarkupList.unique()
			
			genesRegionMap[nomen] = nomenRegion
		}
	} // each line

	return [markNomenMap, sizeMarkMap, genesRegionMap, nomenMarkupMap]
} // loadFeatureNomenclatureFile

/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'binBBFasta.groovy [options] ',
      header:'Options:')
    cli.help('print this message')
    cli.i(longOpt:'in', args:1, argName:'input', 'input uncompressed FASTA file',
      required: true)
    cli.o(longOpt:'out', args:1, argName:'output folder', 'output folder for uncompressed FASTA files',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
