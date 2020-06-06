#!/usr/bin/env groovy

/*
 * annotateMarkup
 * 
 * Annotates the marker string with haplotype sequence features.The 
 * description of the markup must match the descriptions of the fasta.
 * Extracts the approximate gene sequences, erroring on extra on the ends.
 * Disambiguation of ambig markup->gene is done by considering the 
 * previous/left/poroximal gene or the presence of 2DL5.
 * 
 * Input features file. Header row, then lines with 4 tsv columns:
 *  1, markup feature [regex]
 *  2. resolution
 *  3. nomenclature (classic)  [not used]
 *  4. gene nomenclature (classic)
 *  e.g., 
 *  markup	resolution	nomenclature (classic)
 *  AE	gene/allele	3DL3		
 *  BAcL	gene/allele	3DL3		
 *  FGHCIJKL	gene	2DL2L3		
 *  MHCIJKL	gene	2DL1/2DS1/2DS2/2DS2-2DS3
 *  3DL3~FGHCIJKLMHCIRLMHCIJKLSCT	region	cA01	2DL3~2DP1~2DL1~3DP1
 *  WXYKLSCIJKLbRLFZGHCIJKLMHCIJKL~3DL2	region	tB01	2DS2~2DL2~2DL5~2DS3S5~2DP1~2DL1~3DP1
 *  ...
 * 
 * Input markup file. No header; one row per fasta sequence. 
 *  3 columns per row:
 *  1. description of the sequence from fasta
 *  2. markup string
 *  3. csv marker locations using 1-based index from sam file
 *  e.g.,  
 *  NC_000019.10
 *  BCDEFGHCIJKLMHCIRLMHCIJKLSCTUVWXYKLSCIJKLFZZZZZZGHCIJKLSCIJK
 *  241,6750,8644,...
 *
 * Multiple output files will be placed in provided folder. Each will be named the same as the input fasta file. It creates one output files for the markup, and one file per gene; label each seq with hap ID and location 
 * of the gene. (e.g., KP420442.1_3DP1_70186-77539)
 *
 * 
 * annotateMarkup.groovy -i <features definition file> -f <dna fasta file> -m <markup file> -o <output directory>
 * e.g., annotateMarkup.groovy -i ~/doc/kir/snp/15/15_probe_locations_v1_features.txt -f ../NC_000019.10.fasta -m ../NC_000019.10_markerHaps.txt -o output
 * 
 * algorithm
 * Search for a feature in the markup (ie annotate a feature)
 * from largest feature to smallest. Then recursively do the same to 
 * the markup to the left and right. Base case is no markup left or
 * the markup contains no features. Append all annotations or exceptions
 * to make the annotation for the full markup string.
 *
 * @author Dave Roe
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
import java.util.regex.*
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3

// thing that probably won't change per run
err = System.err
//nonGreedyPattern = '^(.*?)'
nonGreedyPattern = '?'
featureSeparator = "~"
annotationExt = "_annotation.txt"
featureFastaExt = "_features.fasta"
// map: gene to approximate UTR sizes
utr5Map = ['3DL3':320, '2DL3':268, '2DP1':267, '3DP1':267, '2DL1':268, '2DL4':267, '3DL1':267, '2DS4':267, '3DL2':268, '2DL2':300, '2DL5':547, '2DS1':268, '2DS2':300, '2DS3':300, '2DS5':268, '3DS1':267, '2DS2-2DS3':300]
utr3Map = ['3DL3':561, '2DL3':510, '2DP1':510, '3DP1':9, '2DL1':510, '2DL4':407, '3DL1':510, '2DS4':624, '3DL2':485, '2DL2':510, '2DL5':415, '2DS1':624, '2DS2':624, '2DS3':645, '2DS5':645, '3DS1':679, '2DS2-2DS3':645]

OptionAccessor options = handleArgs(args)

ArrayList<HashMap> markerRefInfo = loadFeatureNomenclatureFile(new FileReader(options.i))
// map of reference markup to nomenclature features
HashMap<String,String> markerNomenMap = markerRefInfo[0]
// map of length of marker strings to the actual marker strings
HashMap<Integer,ArrayList<String>> sizeMarkerListMap = markerRefInfo[1]
// map of region to gene list; e.g., cA01 -> 3DL3~2DL3~2DP1~2DL1~3DP1
HashMap<String, String> genesRegionMap = markerRefInfo[2]
// map of (gene) nomenclature to Set of markup strings
HashMap<String, TreeSet<String>> nomenMarkupMap = markerRefInfo[3]
// fasta descriptions from the markup sequence; should match the fasta descriptions
ArrayList<String> descList = null
// map of description -> markup
HashMap<String,String> markupMap = null
// map description -> ordered set of location of the probes in the dna
HashMap<String,List<Integer>> locationsMap = null
(descList, descMarkupMap, locationsMap) = loadMarkup(new FileReader(options.m))
if(debugging <= 2) { 
	err.println "locationsMap:" + locationsMap
}

// map description -> dna sequence
//todo: don't store the fasta in memory
HashMap<String,String> faMap = loadFasta(new FileReader(options.f))

if(debugging <= 3) {
	err.println "${faMap.keySet().size()} sequences to markup"
	err.println "${descList.size()} markups with which to interpret"
}

// description -> classic nomenclature based on the markup
HashMap<String,String> descNomenMap = null
// annotate the markup wrt the classic nomenclature
descNomenMap = annotateAll(descMarkupMap, markerNomenMap, sizeMarkerListMap)
if(debugging <= 2) { 
	err.println descNomenMap
}

// populate the description/gene mapping to dna sequences for the gene
// may be more than one sequence per description/gene: CNV
HashBasedTable<String,String,ArrayList<String>> descGeneSeqTable = null
// rows are descriptions, columns are genes, and cells are dna start indexes
HashBasedTable<String,String,ArrayList<String>> descGeneStartIndexTable = null
// rows are descriptions, columns are genes, and cells are dna end indexes
HashBasedTable<String,String,ArrayList<String>> descGeneEndIndexTable = null

(descGeneSeqTable, descGeneStartIndexTable, descGeneEndIndexTable) =
	extractDNAFeaturesFromAll(descMarkupMap, descNomenMap, 
							  nomenMarkupMap, faMap, locationsMap)

if(debugging <= 2) {
	err.println "descriptions and genes in the desc/gene seq table"
	err.println descGeneSeqTable.rowKeySet().join(", ")
	err.println descGeneSeqTable.columnKeySet().join(", ")
}
// output
writeOutput(options.o, options.f, descGeneSeqTable, descGeneStartIndexTable,
			descGeneEndIndexTable, descNomenMap, genesRegionMap)

// end main

/*
 * writeOutput
 */
void writeOutput(String outDir, String fastaName,
				 HashBasedTable<String,String,String> descGeneSeqTable,
				 HashBasedTable<String,String,String> descGeneStartIndexTable,
				 HashBasedTable<String,String,String> descGeneEndIndexTable,
				 HashMap<String,String> descNomenMap, 
				 HashMap<String, String> genesRegionMap) { 
	if(debugging <= 1) {
		err.println "writeOutput(outDir=${outDir}, fastaName=${fastaName})"
	}
	columnSet = descGeneSeqTable.columnKeySet()
	rowSet = descGeneSeqTable.rowKeySet()
	fileName = fastaName.trim().substring(fastaName.lastIndexOf(File.separator)+1, fastaName.length())
	if(debugging <= 2) {
		err.println "writeOutput: fileName=${fileName}"
		//err.println "writeOutput: " + descGeneSeqTable.get("KP420439.1", " 2DL1S1S2")//todo
	}
	//err.println descGeneSeqTable //todo
	// write fastas
	columnSet.each { geneName ->
		//err.println "writeOutput: orig geneName=${geneName}"//todo
		outGeneName = geneName.replaceAll("/2D", "").replaceAll("/3D", "")
		outputFastaName = outDir + "/" +
			fileName.replaceFirst(".fasta", "").replaceFirst(".fa", "") +
			"_" + outGeneName + featureFastaExt
		if(debugging <= 1) {
			err.println "writing ${rowSet.size()} sequences to ${outputFastaName} ..."
		}
		outFastWriter = new PrintWriter(new File(outputFastaName).newOutputStream(), true)
		rowSet.each { desc ->
			//todo: have to make cell a list
			seqs = descGeneSeqTable.get(desc, geneName)
			seqs.eachWithIndex { seq, i ->
				// add locations in the larger string
				start = descGeneStartIndexTable.get(desc, geneName)[i]
				end = descGeneEndIndexTable.get(desc, geneName)[i]
				if(debugging <= 1) {
					err.println "writeOutput: ${desc}, ${geneName} ${start}-${end}"
				}
				//todo: check end - start distance
				outFastWriter.println ">${desc}_${geneName}_${start}-${end}"
				outFastWriter.println seq
			} // each sequence for this gene
		} // each description
		outFastWriter.close()
	} // each gene

	// write annotation
	//todo: make the names match the ones in the feature files
	outputAnnotationName = outDir + "/" +
		fileName.replaceFirst(".fasta", "").replaceFirst(".fa", "") +
		annotationExt
	if(debugging <= 1) {
		err.println "writing to ${outputAnnotationName} ..."
	}
	outAnnWriter = new PrintWriter(new File(outputAnnotationName).newOutputStream(), true)
	rowSet.each { desc ->
		annotation = descNomenMap[desc]
		outAnnWriter.println "${desc}\t${annotation}"
	} // each description
	outAnnWriter.close()
	if(debugging <= 1) {
		err.println "writeOutput: return"
	}
} // writeOutput

/*
 * extractDNAFeaturesFromAll
 * 
 * Given metadata and a nomenclature interpretation of a markup string, 
 * extract the DNA sequences from the fasta.
 * 
 * Given a markup string and its gene annotation, find the location of each
 * gene in the markup and translate the markup positions to DNA positions.
 * Then, extract the DNA of the gene in a generous and inexact manner.
 * 
 * @param descMarkupMap Map of description -> markup
 * @param descNomenMap Map markup -> classic nomenclature
 * @param nomenMarkupMap Map of (gene) nomenclature to Set of markups
 * @param faMap Map description -> dna sequence
 * @param locationsMap HashMap map description -> ordered set of location of the probes in the dna
 * 
 */
ArrayList<HashBasedTable> extractDNAFeaturesFromAll(HashMap<String,String> descMarkupMap,
													HashMap<String,String> descNomenMap,
													HashMap<String,TreeSet<String>> nomenMarkupMap,
													HashMap<String,String> faMap,
													HashMap<String,List<Integer>> locationsMap) {
	if(debugging <= 1){
		err.println "extractDNAFeaturesFromAll()"
	}
	// rows are descriptions, columns are genes, and cells are DNA sequences
	HashBasedTable<String,String,String> descGeneSeqTable = HashBasedTable.create()
	// rows are descriptions, columns are genes, and cells are dna start indexes
	HashBasedTable<String,String,String> descGeneStartIndexTable = HashBasedTable.create()
	// rows are descriptions, columns are genes, and cells are dna end indexes
	HashBasedTable<String,String,String> descGeneEndIndexTable = HashBasedTable.create()

	// each haplotype
	descMarkupMap.each { desc, fullMarkup ->
		workingMarkup = fullMarkup
		workingMarkupIndex = 0 // start from the far left of the markup
		dnaFasta = faMap[desc]
		if(dnaFasta == null) {
			err.println "extractDNAFeaturesFromAll: ERROR: couldn't find DNA sequence for '${desc}'"
			err.println "extractDNAFeaturesFromAll: descriptions from fasta: " + faMap.keySet().join(", ")
			System.exit(1)
		}
		nomenStr = descNomenMap[desc]
		List<Integer> dnaLocationsList = locationsMap[desc]

		if(debugging <= 3) {
			err.println "extractDNAFeaturesFromAll: starting desc=${desc}, fullMarkup=${fullMarkup} (length ${fullMarkup.length()})"
			err.println "fasta size =" + dnaFasta.length()
			//err.println "extractDNAFeaturesFromAll: descNomenMap=${descNomenMap}"
			err.println "extractDNAFeaturesFromAll: dnaLocationsList=" + dnaLocationsList
		}
		if((fullMarkup.length()+1) != dnaLocationsList.size()) {
			err.println "ERROR: markup size doesn't match locations "
			err.println "extractDNAFeaturesFromAll: in desc=${desc}, fullMarkup=${fullMarkup} (length ${fullMarkup.size()})"
			err.println "extractDNAFeaturesFromAll: dnaLocationsList=" + dnaLocationsList
			System.exit(1)
		}
		// todo: remove previousGeneNomen; this is single gene input
		String previousGeneNomen = ""
		nomenStr.split(featureSeparator).each { geneNomen ->
			if(debugging <= 2) {
				err.println "extractDNAFeaturesFromAll: geneNomen=${geneNomen}, workingMarkup=${workingMarkup}"
			}

			// try to find partial 3DL3s and 3DL2s
			// geneNomen can be markup instead of a gene
			//todo if(geneNomen == workingMarkup) {
			if(workingMarkup.startsWith(geneNomen)) {
				// check for partal 3DL3
				geneMarkupSet = nomenMarkupMap["3DL3"]
				if(debugging <= 2) {
					err.println "3DL3 geneMarkupSet=${geneMarkupSet}"
				}

				found = false
				for(Iterator gmIter = geneMarkupSet.iterator();
					(gmIter.hasNext() && !found);) {
					geneMarkup = gmIter.next()
					//todo if(geneMarkup.endsWith(workingMarkup)) {
					if(geneMarkup.endsWith(geneNomen)) { 
						if(debugging <= 3) {
							err.println "extractDNAFeaturesFromAll: found 3DL3(partial) ${workingMarkup} in ${geneMarkup}"
						}
						(workingMarkupNew, workingMarkupIndexNew) =
							extractDNA(desc, previousGeneNomen, nomenStr,
									   geneNomen, "3DL3", true,
									   workingMarkup, workingMarkupIndex,
									   dnaLocationsList, dnaFasta, descGeneSeqTable,
									   descGeneStartIndexTable, descGeneEndIndexTable)
						if(debugging <= 2) { 
							err.println "${workingMarkupNew}, ${workingMarkupIndexNew}"
						}
						workingMarkup = workingMarkupNew
						workingMarkupIndex = workingMarkupIndexNew
						previousGeneNomen = geneNomen
						found = true
					} // found 3DL3
				} // each 3DL3 allele
				// check for partal 3DL2
				geneMarkupSet = nomenMarkupMap["3DL2"]
				if(debugging <= 2) {
					err.println "3DL2 geneMarkupSet=${geneMarkupSet}"
				}
				found = false
				for(Iterator gmIter = geneMarkupSet.iterator();
					(gmIter.hasNext() && !found);) {
					geneMarkup = gmIter.next()
					//todo(old) if(geneMarkup.startsWith(workingMarkup)) {
					if(geneMarkup.startsWith(geneNomen)) { 
						if(debugging <= 3) {
							err.println "extractDNAFeaturesFromAll: found 3DL2(partial) ${workingMarkup} in ${geneMarkup}"
						}
						(workingMarkupNew, workingMarkupIndexNew) =
							extractDNA(desc, previousGeneNomen, nomenStr,
									   geneNomen, "3DL2", true,
									   workingMarkup, workingMarkupIndex,
									   dnaLocationsList, dnaFasta, descGeneSeqTable,
									   descGeneStartIndexTable, descGeneEndIndexTable)
						if(debugging <= 2) {
							err.println "${workingMarkupNew}, ${workingMarkupIndexNew}"
						}
						workingMarkup = workingMarkupNew
						workingMarkupIndex = workingMarkupIndexNew
						previousGeneNomen = geneNomen
						found = true
					} // found 3DL2
				} // each 3DL2 allele
				return // next geneNomen
			} // doesn't match a gene; potential partial 3DL3 or 3DL2
			
			// loop through all the markers for each geneNomen
			geneMarkupSet = nomenMarkupMap[geneNomen]

			for(Iterator gmIter = geneMarkupSet.iterator(); gmIter.hasNext();) {
				geneMarkup = gmIter.next()
				(workingMarkupNew, workingMarkupIndexNew) =
					extractDNA(desc, previousGeneNomen, nomenStr, geneMarkup,
							   geneNomen, false,
							   workingMarkup, workingMarkupIndex,
							   dnaLocationsList, dnaFasta, descGeneSeqTable,
							   descGeneStartIndexTable, descGeneEndIndexTable)
				if(workingMarkupNew != "") {
					if(debugging <= 3) {
						err.println "extractDNAFeaturesFromAll: found ${geneNomen} ${geneMarkup} in ${workingMarkup}"
						err.println "extractDNAFeaturesFromAll: workingMarkupIndex=${workingMarkupIndex} workingMarkupIndexNew=${workingMarkupIndexNew}"
					}
					workingMarkup = workingMarkupNew
					workingMarkupIndex = workingMarkupIndexNew
					previousGeneNomen = geneNomen
					break // break out of the markup for this gene
				}
			} // each markup for this gene
			if(workingMarkup == "") {
				err.println "extractDNAFeaturesFromAll: ERROR: couldn't find ${geneNomen} markup in part of full markup ${fullMarkup}"
				System.exit(1)
			}
		} // each genNomen
	} // each description and its markup

	if(debugging <= 1) {
		err.println "extractDNAFeaturesFromAll: return ${descGeneSeqTable.rowKeySet().size()} descriptions and ${descGeneSeqTable.columnKeySet().size()} genes"
	}
	return [descGeneSeqTable, descGeneStartIndexTable, descGeneEndIndexTable]
} // extractDNAFeaturesFromAll

/*
 * extractDNA
 *
 * @param desc String of the name of the sequence being annotated
 * @return the new workingMarkup and workingMarkupIndex or "" and 0 on error
 */
ArrayList extractDNA(String desc, String previousGeneNomen, String nomenStr,
					 String geneMarkup,
					 String geneNomen, Boolean partial, String workingMarkup,
					 Integer workingMarkupIndex,
					 ArrayList<Integer> dnaLocationsList,
					 String dnaFasta,
					 HashBasedTable<String,String,String> descGeneSeqTable,
					 HashBasedTable<String,String,String> descGeneStartIndexTable,
					 HashBasedTable<String,String,String> descGeneEndIndexTable) {
	if(debugging <= 1) {
		err.println "extractDNA(desc=${desc}, previousGeneNomen=${previousGeneNomen}, nomenStr=${nomenStr}, geneMarkup=${geneMarkup}, geneNomen=${geneNomen}, partial=${partial}, workingMarkup=${workingMarkup}, workingMarkupIndex=${workingMarkupIndex})"
	}
	// use a regex to find the markup feature in the working markup
	Pattern pattern = Pattern.compile(geneMarkup)
	Matcher matcher = pattern.matcher(workingMarkup)
	boolean found = matcher.find()
	if(debugging <= 2) { 
		err.println "extractDNA: found=${found}"
	}
	if(found == false) {
		if(debugging <= 1) {
			err.println "extractDNA: couldn't find ${geneMarkup} in ${workingMarkup}"
			err.println "extractDNA: return ['', 0]"
		}
		return ["", 0]
	}
	/* e.g. 
	 * BCDEFGHCIJKLMHCIRLMHCIJKLSCTUVWXYKLSCIJKLFZZZZZZGHCIJKLSCIJK
	 * 241,6750,8644,12187,17578,18063,18892,21119,22947,26401,27658,31713,34055,35499,37747,39583,42957,46541,48883,50234,52472,54316,57678,58935,62983,65387,67337,69176,72740,76339,81857,85102,88323,89159,93184,95533,98852,100716,104094,105350,109445,111797,112317,112451,112682,112778,112913,113009,113067,114464,116715,118579,121934,123191,127248,129508,132886,134770,138145,139402,145868
	 */
	Integer markupIndex = matcher.start()
	Integer markupEndIndex = matcher.end()
	int i = markupIndex + workingMarkupIndex
	int ie = markupEndIndex + workingMarkupIndex
	if(debugging <= 2) {
		err.println "extractDNA: ${dnaLocationsList.size()} items in ${dnaLocationsList}"
		err.println "extractDNA: markupIndex=${markupIndex}"
		err.println "extractDNA: markupEndIndex=${markupEndIndex}"
		err.println "extractDNA: workingMarkupIndex=${workingMarkupIndex} (${dnaLocationsList[markupEndIndex]})"
		err.println "i=${i}, ie=${ie}"
		err.println "extractDNA: dnaLocationsList[{$i}]=${dnaLocationsList[i]}"
		err.println "extractDNA: dnaLocationsList[${markupEndIndex} + ${workingMarkupIndex}]=" + dnaLocationsList[ie]
	}
	Integer dnaIndex = dnaLocationsList[i]
	Integer dnaEndIndex = dnaLocationsList[ie]
	if(partial == true) {
		if(geneNomen == "3DL3") {
			dnaIndex = 0
		} else if(geneNomen == "3DL1L2"){
//            dnaIndex -= 5000 // not for MN167507
			dnaEndIndex = dnaFasta.length() - 1
		}
	} else if(geneNomen == "2DL5") {
        // needs more padding on proximal(centromeric) end
        dnaIndex -= 3500
    } else if(geneNomen == "3DL1S1") {
        // needs more padding on proximal(centromeric) end
        dnaEndIndex += 2500
    } else if(geneNomen == "3DL1L2") { // for the fusion
        // needs more padding on proximal(centromeric) end
        // and distal (telomeric) end
//todo(put bak? test wih MN167530)        dnaIndex -= 200
        dnaIndex -= 2400
        dnaEndIndex += 2000
    } else if(geneNomen == "3DL2") {
        // needs more padding on proximal(centromeric) end
        // and distal (telomeric) end
//todo(put bak? test wih MN167530)        dnaIndex -= 200
        dnaIndex -= 5000
        dnaEndIndex += dnaFasta.length() - 1
    } else if(geneNomen == "3DP1") {
        // needs more padding on distal(telomeric) end
        dnaEndIndex += 2000
    } else if((geneNomen == "2DL2L3S3S4S5") || (geneNomen == "2DS4") 
              || (geneNomen == "2DL2")) {
//    } else if((geneNomen == "2DL2L3S3S4S5")) {
        // 2DS4 needs more padding on both ends
        dnaIndex -= 2000 // MN167521 2DS4 and MN167525 2DL2
        dnaEndIndex += 1500 // sometimes shouldn't be used; sometimes should (MN16752521)
    }
	if(debugging <= 2) {
		err.println "extractDNA: dnaIndex=${dnaIndex}"
		err.println "extractDNA: dnaEndIndex=${dnaEndIndex}"
	}
	dnaEndIndexStore = dnaEndIndex
	if((dnaEndIndex == -1) || (dnaEndIndex > dnaFasta.length()-1)) {
		dnaEndIndexStore = dnaFasta.length() - 1
        dnaEndIndex = dnaFasta.length() - 1
	}
	geneDNA = dnaFasta[dnaIndex..dnaEndIndex]
	haplotypeLocus = geneToHaplotypeLocus(previousGeneNomen, nomenStr,
										  geneNomen)
	seqList = descGeneSeqTable.get(desc, haplotypeLocus)
	if(seqList == null) {
		seqList = new ArrayList(1)
		startList = new ArrayList(1)
		endList = new ArrayList(1)
		seqList.add(geneDNA)
		startList.add(dnaIndex)
		endList.add(dnaEndIndexStore)
		descGeneSeqTable.put(desc, haplotypeLocus, seqList)
		descGeneStartIndexTable.put(desc, haplotypeLocus, startList)
		descGeneEndIndexTable.put(desc, haplotypeLocus, endList)
	}
    else {
		if(debugging <= 4) {
			err.println "extractDNA: multiple sequences for ${haplotypeLocus} in ${desc}"
//            err.println "keeping the longer one"// todo: maybe not the best?
		}
  // this doesn't work with multiple different genes (e.g., 2DL2L3S3S4S5)
  // but is needed for some genes with g1 g2 etc from augustus
  /*        String firstSeq = seqList[0]
          if(firstSeq.length() < (dnaEndIndexStore - dnaIndex)) {
		    seqList[0] = geneDNA
            descGeneStartIndexTable.get(desc, haplotypeLocus).clear()
            descGeneEndIndexTable.get(desc, haplotypeLocus).clear()
*/
		seqList.add(geneDNA)
		startList.add(dnaIndex)
		endList.add(dnaEndIndexStore)
		descGeneSeqTable.put(desc, haplotypeLocus, seqList)
		    descGeneStartIndexTable.get(desc, haplotypeLocus).add(dnaIndex)
		    descGeneEndIndexTable.get(desc, haplotypeLocus).add(dnaEndIndexStore)
//        }
	}

	if(debugging <= 2) {
		err.println "extractDNA: putting value in descGeneSeqTable: ${desc}, ${haplotypeLocus} = DNA of size ${geneDNA.length()}"
		err.println "extractDNA: ${dnaIndex} - ${dnaEndIndexStore}"
		err.println "extractDNA: start index(${desc}, ${haplotypeLocus})=" +
			descGeneStartIndexTable.get(desc, haplotypeLocus)
		err.println "extractDNA: end index(${desc}, ${haplotypeLocus})=" +
			descGeneEndIndexTable.get(desc, haplotypeLocus)
	}
	Integer length = dnaEndIndexStore - dnaIndex
	if(length > 30000) {
		err.println "extractDNA: WARNING: suspicious length (${length}) for ${desc}, ${haplotypeLocus}"
		//System.exit(1)//todo (remove?)
	}
	// now work on the rest of the markup to the right
	workingMarkupIndex = ie
	if(debugging <= 2) {
		err.println "extractDNA: markupEndIndex=${markupEndIndex}, workingMarkupIndex=${workingMarkupIndex}, workingMarkup length=${workingMarkup.length()}"
	}
	if(markupEndIndex < workingMarkup.length()) {
		workingMarkup = workingMarkup[markupEndIndex..-1]
	} else {
		workingMarkup = ""
		workingMarkupIndex = 0
	}

 	if(debugging <= 1) {
		err.println "extractDNA: return [${workingMarkup}, ${workingMarkupIndex}]"
	}
	return [workingMarkup, workingMarkupIndex]
} // extractDNA


/*
 * geneToHaplotypeLocus
 * 
 * Converts the ambiguous markup gene names to haplotype locus (IPD-KIR) names.
 * This is done by considering the previous/left/poroximal gene.
 * 2DL1/2DS1/2DS2/2DS2\2DS3: based on 3DL3, 2DS, and 2DP1
 * 2DS3S5/2DS4, 2DL2L3, 2DL1S1: based on 2DL5
 */
String geneToHaplotypeLocus(String previousGeneNomen, String nomenStr,
							String geneNomen) {
	if(debugging <= 1) {
		err.println "geneToHaplotypeLocus(previousGeneNomen=${previousGeneNomen}, geneNomen=${geneNomen})"
	}
	String ret = geneNomen
	if(geneNomen.contains("2DL1")) { // 2DL1/2DS1/2DS2/2DS2\2DS3
        ret = "2DL1S1S2"
	} else if(geneNomen.contains("2DS3")) { // 2DS3S5/2DS4
        ret = "2DS3S4S5"
	} /*else if(geneNomen.contains("2DL2L3")) { // 2DL2L3
        if(previousGeneNomen.contains("2DS2")) {
            ret = "2DL2"
		} else if(!previousGeneNomen.contains("2DS2")) {
			ret = "2DL3"
		} else {
			err.println "geneToHaplotypeLocus: WARNING: couldn't find previous gene for previousGeneNomen=${previousGeneNomen}, nomenStr=${nomenStr}, geneNomen=${geneNomen})"
			ret = "2DL3"
		}
	} else if(geneNomen.contains("3DL1S1")) { // 3DL1S1
        if(nomenStr.contains("2DL5")) {
            ret = "3DS1"
		} else if(!nomenStr.contains("2DL5")) {
			ret = "3DL1"
		} else {
			err.println "geneToHaplotypeLocus: WARNING: couldn't find previous gene for previousGeneNomen=${previousGeneNomen}, nomenStr=${nomenStr}, geneNomen=${geneNomen})"
		}
	}*/
	if(debugging <= 1) {
		err.println "geneToHaplotypeLocus: return ${ret}"
	}
	return ret
} //geneToHaplotypeLocus

/*
 * annotate
 * 
 * Annotates a markup strings with KIR haplotypes and gene features.
 * 
 * @param markupMap Map of description -> markup (string representing the probes along the dna string)
 * @param markerNomenMap a Map of reference marker strings to KIR features
 * @param sizeMarkerListMap a Map of sizes of reference features to Lists of features strings of that size
 * @return Map of descriptions to each's nomenclature call
 */
HashMap<String,String> annotateAll(HashMap<String,String> markupMap,
								   HashMap<String,String> markerNomenMap,
								   HashMap<Integer,ArrayList<String>> sizeMarkerListMap) {
	if(debugging <= 1){
		err.println "annotateAll(${markupMap})"
	}
	// description -> classic nomenclature based on the markup
	HashMap<String,String> descNomenMap = new HashMap()
	
	for(Iterator muIter = markupMap.keySet().iterator(); muIter.hasNext();) {
		desc = muIter.next()
		markup = markupMap[desc]
		annotation = annotate(markup, markerNomenMap, sizeMarkerListMap)
		descNomenMap[desc] = annotation
	} // each markup string
		
	if(debugging <= 1){
		err.println "annotateAll: return"
	}
	return descNomenMap
} // annotateAll

/*
 * Annotate the markup with the info from information about the markers.
 * Populate a table that stores the DNA sequences for each descriiption
 * and gene. Also, a map that stores the (classic) nomenclature 
 * for the description.
 *
 * @return the (classic) nomenclature annotation of the markup
 */
// what to return?
String annotate(String markup,
				HashMap<String,String> markerNomenMap,
				HashMap<Integer,ArrayList<String>> sizeMarkerListMap) {
    if(debugging <= 1) {
		err.println "annotate(markup=${markup})"
	}
	if((markup == null) || (markup == "")) { // base case
        if(debugging <= 1){
			err.println "annotate: return \"\""
		}
		return ""
	}
	// iterate through the markup features, from longest to shortest
	List<Integer> sizeList = new ArrayList<Integer>(sizeMarkerListMap.keySet())
	sizeList = sizeList.sort().reverse()
	if(debugging <= 1) { 
		err.println "sizeList = " + sizeList
	}
	for(Iterator siter = sizeList.iterator(); siter.hasNext();) {
		Integer size = siter.next()
		// iterate through each markup feature of this size
		featureList = sizeMarkerListMap[size]
		for(Iterator fiter = featureList.iterator(); fiter.hasNext();) {
			String feature = fiter.next()
			String nomen = markerNomenMap[feature]
			if(debugging <= 1) {
				err.println "annotate: markup=${markup}"
				err.println "annotate: checking ${feature} (${nomen})"
			}

			// use a regex to find the feature
			//todo
			Pattern pattern = Pattern.compile(feature)
//todo			Pattern pattern = Pattern.compile("(${feature})${nonGreedyPattern}")
			Matcher matcher = pattern.matcher(markup)
			boolean found = matcher.find()
			//todo boolean found = matcher.matches()
			//err.println "annotate: find $feature: ${found}"//todo
			if(found == false) { //can't trust this; see continue below
				continue;
			}
			//err.println "droe=${matcher.start()}"//droe
			Integer fIndex = matcher.start()
			Integer fendIndex = matcher.end()

			if(fIndex == fendIndex) { // didn't find the pattern
				continue;
			}
			// todo: markupPre isn't getting set correct
			String markupPre = ""
			if(fIndex > 0) { 
				markupPre = markup[0..fIndex-1]
			}

			int postStart = fIndex + (fendIndex - fIndex)//todo
			String markupPost = ""
			if(debugging <= 1) {
				err.println "markupPre=$markupPre"
				err.println "markup length=${markup.length()}"
				err.println "postStart=$postStart, length=${markup.length()}"
			}
			if(postStart < markup.length()) {
				markupPost = markup[postStart..-1]
			}

			if(debugging <= 1){
				err.println "annotate: found ${feature}, fIndex=${fIndex}"
				err.println "feature=${feature}, nomenclature=${nomen}"
				err.println "annotate: markupPre=${markupPre}, markupPost=${markupPost}"
			}
			// recursively call and return
			String preResult = annotate(markupPre, markerNomenMap,
										   sizeMarkerListMap)
			String postResult = annotate(markupPost, markerNomenMap,
											sizeMarkerListMap)
			if(debugging <= 3) {
				err.println "preResult=${preResult} (markup=${markup})"
				err.println "postResult=${postResult} (markup=${markup})"
				err.println "nomenclature=${nomen} (markup=${markup})"
			}
			if(preResult != "") {
				preResult = preResult + featureSeparator
			}
			if(postResult != "") {
				postResult = featureSeparator + postResult
			}
			ret = preResult + nomen + postResult
			if(debugging <= 1) {
				err.println "annotate(${markup}): return ${ret}"
			}
			return preResult + nomen + postResult
			
		} // each feature of this size
	} // each size, from

	// the markup is not a feature, return the markup; base case
	if(debugging <= 1){
		err.println "annotate: no nomenclature; return ${markup}"
	}
	return markup
} // annotate

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
ArrayList<HashMap> loadFeatureNomenclatureFile(FileReader reader) {
	// markup string -> nomenclature string
    HashMap<String,String> markNomenMap = new HashMap()
	// integer size -> list of markup strings
    HashMap<Integer,ArrayList<String>> sizeMarkMap = new HashMap()
	// map of region to gene list; e.g., cA01 -> 3DL3~2DL3~2DP1~2DL1~3DP1
	HashMap<String, String> genesRegionMap = new HashMap()
	// map of (gene) nomenclature to markup
	HashMap<String,TreeSet<String>> nomenMarkupMap = new HashMap()
	
	reader.readLine() // header
	reader.eachLine { line ->
		if(debugging <= 1) {
			err.println "loadFeatureNomenclatureFile: line=${line}"
		}
        ArrayList cols = line.split('\t')
        String mark = cols[0].trim() // markup
		//mark = mark.replaceFirst("3DL3~", "").replaceFirst("~3DL2", "")
        String resolution = cols[1].trim() // resolution
		String nomenRegion = cols[2].trim()
        String nomen = cols[3] ? cols[3].trim() : null // gene nomanclature (classic)
		if(nomen == null) {
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
		nomenMarkupMap[nomen] = nomenMarkupList
		
		if(resolution == "region") {
			genesRegionMap[nomen] = nomenRegion
		}
/*		if(!nomen.contains("2DL3")) { // testing
			return
		}
*/
	} // each line

	return [markNomenMap, sizeMarkMap, genesRegionMap, nomenMarkupMap]
} // loadFeatureNomenclatureFile

/*
 * loadMarkup
 * Load the markup file and return its descriptions, markups, and locations.
 * See the top-level documentation for the input format.
 *
 * @return three Collections
 *  ArrayList of descriptions of the markup sequence
 *  Map of description -> markup
 *  Map of description -> ordered set of location of the probes in the dna
 *
 */
List loadMarkup(FileReader reader) {
	// descriptions of the markup sequence; should match the fasta descriptions
	ArrayList<String> descList = new ArrayList()
	// map of description -> markup
	HashMap<String,String> descMarkupMap = new HashMap()
	// map description -> ordered set of location of the probes in the dna
	HashMap<String,List<Integer>> locationsMap = new HashMap()

	reader.eachLine { line ->
		if(debugging <= 1) {
			err.println "loadFeatureNomenclatureFile: line=${line}"
		}
		// 'probe pairs' and things in brackets are informational
		if((line == null) || (line == "") ||
		   line.contains("probe pairs") || line.contains("[")) {
			return
		}
		ArrayList cols = line.split('\t')
		String description = cols[0]
		String markup = cols[1]
		String locationString = cols[2]
		descList.add(description)
		descMarkupMap[description] = markup
		// convert csv to 
		ArrayList<Integer> locationList =
			locationString.split(',').collect { it as int}
		locationsMap[description] = locationList
	} // each lin

	return [descList, descMarkupMap, locationsMap]
} // loadMarkup

/*
 * loadFasta
 * Load the fasta file and return a map of descriptions to DNA strings.
 * 
 * @param FileReader fasta file of DNA strings
 * @return Map of the description -> DNA string
 *
 * e.g., 
 *  >NC_000019.10
 *  TCCTAAGTGAACTAACC...
 */
HashMap<String,String> loadFasta(FileReader reader) {
    HashMap<String,String> faList = new HashMap()
	String desc = null
	StringBuilder dna = new StringBuilder()

	reader.eachLine { line ->
		if(line.startsWith(">")) {
			if(desc != null) {
				// set the old value
				faList[desc] = dna.toString()
				// reset for the current line
				desc = null
				dna = new StringBuilder()
			}
			desc = line.trim().replaceFirst(">", "")
			// truncate to space
			// todo: paramaterize this; coordinate with alignment2ProbePairs
			descArray = desc.split(' ')
			desc = descArray[0]
			return
		}
		if(desc == null) {
			return
		}
		dna.append(line.trim())
	} // each line
	if(desc != null) {
		faList[desc] = dna.toString()
	}

	if(debugging <= 1) {
		err.println "loadFasta: return=${faList.size()} fasta sequences"
	}
	return faList
} // loadFasta

/*
 * handleArgs

 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'alignment2ProbePairs.groovy [options] ', header:'Options:')
    cli.i(longOpt:'input', args:1, argName:'in', 'input of features and their nomenclature',
		  required: true)
    cli.f(longOpt:'fasta', args:1, argName:'fasta', 'fasta containing the haplotype DNA sequence',
		  required: false)
    cli.m(longOpt:'markup', args:1, argName:'markup', 'text file containing the probe markup string for the fasta',
		  required: false)
    cli.o(longOpt:'directory to put the output', args:1, argName:'out', 
		  'output directory', required: true)
    OptionAccessor options = cli.parse(args)

    return options
} // handleArgs
