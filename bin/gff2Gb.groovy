#!/usr/bin/env groovy

/*
 * gff2Gb
 * 
 * Converts the output from Gaius-Augustus to input suitable for NCBI's
 * table2asn_GFF.
 * 
 * https://github.com/Gaius-Augustus/Augustus
 * https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
 * 
 * gff2Gb.groovy -i todo
 * e.g., gff2Gb.groovy -d KP420440 -i KP420440_2DL2L3S3S4S5_augustus.gff -l KP420440_2DL2L3S3S4S5.gl.txt -o KP420440_2DL2L3S3S4S5.gff 2> KP420440_2DL2L3S3S4S5_gffGb_err.txt
 * 
 * @author Dave Roe
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
import groovy.transform.Field

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3

// thing that probably won't change per run
err = System.err
// type values in gff
@Field final HashMap<String, String> childTypeMap = ['gene':'mRNA', 'mRNA':'CDS']
// region value to ID abbreviations
@Field final HashMap<String, String> regionAbbrevMap = ['mRNA':'mRNA', 'five_prime_UTR':'5utr', 'CDS':'CDS', 'three_prime_UTR':'3utr']
/*
 * from the GFF spec (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) 
 *
 * Column 9: "attributes"
 *    A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons.
 * ID, Name, Alias, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term, Is_circular
 */
@Field final Integer SEQID_COL_INDEX = 0
@Field final Integer TYPE_COL_INDEX = 2
@Field final Integer START_COL_INDEX = 3
@Field final Integer END_COL_INDEX = 4
@Field final Integer ATTRIBUTE_COL_INDEX = 8

@Field final OptionAccessor options = handleArgs(args)

FileReader iReader = new FileReader(options.i)
FileReader lReader = new FileReader(options.l)
@Field final PrintWriter writer = new PrintWriter(new File(options.o).newOutputStream(), true)
if(debugging <= 3) {
	err.println "reading from ${options.i} and ${options.l}, writing to ${options.o} ..."
}

// get the locus
iArray = options.i.split('_')
locus = "KIR" + iArray[--2]

// write the gff meta data
writer.println "##gff-version 3"
writer.println "#!gff-spec-version 1.21"

// Map of sequence descriptions to their AUGUSTUS ID (e.g., g1.t1)
HashMap<String,String> seqAugNameMap = new HashMap<String,String>()
// Map of sequence descriptions to their best IPD-KIR allele name (e.g., 0030101)
HashMap<String,String> seqGLMap = new HashMap<String,String>()
//  Map of sequence descriptions to resolution of best IPD-KIR alleles name
HashMap<String,String> seqGLResMap = new HashMap<String,String>()
loadGL(lReader, seqAugNameMap, seqGLMap, seqGLResMap)

processType("gene", -1, "", iReader, seqAugNameMap, seqGLMap, seqGLResMap,
            options.d, locus)

if(debugging <= 3) {
	err.println "done"
}
iReader.close()
lReader.close()
writer.close()

/*
 * processType
 *
 * @params inType String indicating the current type being processed
 * @params inTypeIndex Integer indicating the index of the type currently being processed; -1 initially
 * @params parentName String the name (including index) of parent (e.g., gene0)
 * @params iReader FileReader for the gff file
 * @params seqAugNameMap Map of sequence descriptions to their AUGUSTUS ID (e.g., g1.t1)
 * @params seqGLMap Map of sequence descriptions to their best IPD-KIR allele name (e.g., 0030101)
 * @params seqGLResMap Map of sequence descriptions to resolution of best IPD-KIR alleles name
 * @param origDesc String the description from the original contig/haplotype fasta (not the feature)
 * @param locus String the locus of the features
 * 
 * Process a gene type and all its children.
 */
void processType(String inType, Integer inTypeIndex, String parentName,
                 FileReader iReader, HashMap<String,String> seqAugNameMap,
                 HashMap<String,String> seqGLMap, HashMap<String,String> seqGLResMap,
                 String origDesc, String locus) {
	if(debugging <= 1) {
		err.println "processType(inType=${inType}, inTypeIndex=${inTypeIndex}, parentName=${parentName}, origDesc=${origDesc}, locus=${locus})"
	}

    Integer typeIndex = inTypeIndex
    String nextType = ""
    while((line = iReader.readLine()) != null) {
        if(debugging <= 1) {
		    err.println "processType: line=${line}"
            err.println "processType: parentName=${parentName} inType=${inType}, typeIndex=${typeIndex}"
	    }

        ArrayList<String> cols = line.split('\t')
        String type = cols[TYPE_COL_INDEX]
        if(isChild(inType, type) == true) { // e.g., mRNA and gene
            if(debugging <= 1) {
		        err.println "processType(flip child): (type=${type}, inType=${inType})"
	        }
            processAttribute(type, typeIndex, cols, writer, parentName, seqAugNameMap,
                             origDesc, locus)
            if(debugging <= 1) {
		        err.println "processType(flip child): return"
	        }
           return
        } else if(type == inType) {
            if(type == "gene") {
                parentName = ""
            }
            err.println "typeIndex2=${typeIndex}2"
            typeIndex++
            err.println "typeIndex3=${typeIndex}"
            processAttribute(type, typeIndex, cols, writer, parentName, seqAugNameMap,
                             origDesc, locus)
//remove(todo)        } else if(type == nextType){
        } else if(isChild(type, inType) == true) {
            Integer parentIndex = typeIndex
            parentName = parentName + inType + parentIndex   // e.g., gene0
            processAttribute(type, 0, cols, writer, parentName, seqAugNameMap,
                             origDesc, locus)
//remove(todo)            if(childTypeMap[type] != null) { // if a child is expected
                processType(type, 0, parentName, iReader, seqAugNameMap, seqGLMap,
                            seqGLResMap, origDesc, locus)
//remove(todo)            }
	        if(debugging <= 2) {
		        err.println "processType(child): type=${type}, parentName=${parentName}"
	        }
            if(childTypeMap[type] == null) {
	            if(debugging <= 1) {
		            err.println "processType: return"
	            }
                return
            }
        } else {
            parentNameUTR = parentName
            if(type == "transcription_start_site") {
                type = "five_prime_UTR"
            } else if (type == "transcription_end_site") {
                type = "three_prime_UTR"
//                (parentNameUTR, rest) = parentName.split("mRNA")
            }
            if(type.contains("UTR")){
                cols[TYPE_COL_INDEX] = type
                processAttribute(type, 0, cols, writer,
                                 parentNameUTR, seqAugNameMap, origDesc, locus)
            }
            if (type == "three_prime_UTR") {
                type = origDesc
                parentName = ""
	            if(debugging <= 1) {
		            err.println "processType(tts): return"
	            }
                return
            }

        }
    } // each line in GFF

	if(debugging <= 1) {
		err.println "processType: return"
	}
} // processType

/*
 * processAttribute
 *
 * Converts the attribute line and writes the line to the output.
 * 
 * @params inType String indicating the current type being processed
 * @params typeIndex Integer indicating the index of the type currently being processed
 * @params cols ArrayList of columns of the gff line that starts this gene
 * @param origDesc String the description from the original contig/haplotype fasta (not the feature)
 * @return String the next (child) type or "" if none
 */
String processAttribute(String inType, Integer typeIndex, ArrayList<String> cols,
                        PrintWriter writer, String parentName,
                        HashMap<String,String> seqAugNameMap, String origDesc,
                        String locus) {
	if(debugging <= 1) {
		err.println "processAttribute(inType=${inType}, typeIndex=${typeIndex}, parentName=${parentName})"
	}

    // replace the first gene ID with the locus
/*    if(parentName == "") {
        parentName = locus
    } else {*/
        parentName = parentName.replaceFirst("gene", locus)
//    }
    
	if(debugging <= 1) {
		err.println "processAttribute: new parentName=${parentName}"
	}
    
    // lift the feature coordinates to the full contig/haplotype sequence
    desc = cols[SEQID_COL_INDEX] // e.g., KP420440_2DL2L3S3S4S5_18152-34641
    coordStartIndex = desc.lastIndexOf('_')
    (coordStart, coordEnd) = desc[coordStartIndex+1..-1].split('-')
    Integer coordStart = coordStart.toInteger()
    Integer coordEnd = coordEnd.toInteger()
	if(debugging <= 2) {
		err.println "processAttribute: desc=${desc}, coordStart=${coordStart}, coordEnd=${coordEnd}"
		err.println "processAttribute: initial feature start=${cols[START_COL_INDEX]}, initial feature end=${cols[END_COL_INDEX]}"
	}
    cols[START_COL_INDEX] = cols[START_COL_INDEX].toInteger() + coordStart
    cols[END_COL_INDEX] = cols[END_COL_INDEX].toInteger() + coordEnd
    String att = cols[ATTRIBUTE_COL_INDEX]
	if(debugging <= 2) {
		err.println "processAttribute: new feature start=${cols[START_COL_INDEX]}, new feature end=${cols[END_COL_INDEX]}"
        err.println "processAttribute: att=${att}"
	}

    // filter the augustus gff by the best gene/transcript
    // as determined in the gl file
    String printType =  inType
    if(printType == "gene") {
        printType = locus
    }
    ArrayList attCols = new ArrayList()
    attCols.addAll(att.split(';'))
    id = attCols[0]
    (rest, currentSeqName) = id.split('=')
    currentSeqName = currentSeqName.replaceFirst("tss", "").replaceFirst("tts", "")
    String bestAugName = seqAugNameMap[desc]
	if(debugging <= 2) {
		err.println "processAttribute: current name from augustus=${currentSeqName}"
        err.println "processAttribute: best gene name from augustus=${bestAugName}"
	}

    String newAtt = "ID=${parentName}${printType}${typeIndex};"
    if(parentName != "") {
        newAtt += "Parent=${parentName}"
    }
    // e.g., g1.t1 and
    // todo: have to change the desc to the full contig
    if(bestAugName.startsWith(currentSeqName) ||
       currentSeqName.startsWith(bestAugName)) {
        outLine = "${origDesc}\t"
        outLine += cols[1..(ATTRIBUTE_COL_INDEX-1)].join('\t')
        outLine += "\t${newAtt}"
	    if(debugging <= 3) {
		    err.println "processAttribute: writing ${outLine}"
        }
        writer.println outLine
    }

    String nextType = childTypeMap[inType] ? childTypeMap[inType] : ""
	if(debugging <= 1) {
		err.println "processAttribute: return ${nextType}"
	}
    return nextType
} // processAttribute

/*
 * isChild
 *
 * Tests if cType is a child of pType.
 */
boolean isChild(String cType, String pType) {
    if(debugging <= 1) {
		err.println "isChild(cType=${cType}, pType=${pType})"
	}
    boolean ret = false
    if((pType == "gene") && (cType == "mRNA")) {
        ret = true
    } else if((pType == "mRNA") && (cType == "CDS")) {
        ret = true
    }
    if(debugging <= 1) {
		err.println "isChild: return ${ret}"
	}
    return ret
} // isChild

/*
 * loadGL
 *
 */
void loadGL(FileReader lReader, HashMap<String,String> seqAugNameMap,
            HashMap<String,String> seqGLMap, HashMap<String,String> seqGLResMap) {
    lReader.eachLine { line ->
        if(line.trim() == "") {
            return
        }
        // ID	augID	gene	protein	cDNA	full
        (id, augID, gene, protein, cDNA, full) = line.split('\t')
        seqAugNameMap[id] = augID
        if(!full.contains("NEW")) {
            bestAllele = full
            bestRes = "full gene"
        } else if(!cDNA.contains("NEW")) {
            bestAllele = cDNA
            bestRes = "cDNA"
        } else if(!protein.contains("NEW")) {
            bestAllele = protein
            bestRes = "protein"
        }
        seqGLMap[id] = bestAllele
        seqGLResMap[id] = bestRes
    }

} // loadGL

/*
 * handleArgs

 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'gff2Gb.groovy [options] ', header:'Options:')
    cli.d(longOpt:'description', args:1, argName:'desc', 'description in the original fasta',
		  required: true)
    cli.i(longOpt:'input', args:1, argName:'in', 'GFF output from AUGUSTUS',
		  required: true)
    cli.l(longOpt:'gl', args:1, argName:'gl', 'gl.txt file',
		  required: true)
    cli.o(longOpt:'modified GFF', args:1, argName:'out', 
		  'output directory', required: true)
    OptionAccessor options = cli.parse(args)

    return options
} // handleArgs
