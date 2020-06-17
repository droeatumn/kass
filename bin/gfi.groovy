#!/usr/bin/env groovy
/*
 * Gene Feature Interpretation
 *
 * Annotates a gene sequences with their UTRs, introns, and exons 
 * and assigns IPD-KIR names. Assumes features cannot be more than maxFeatureDistance bp away.
 *
 * e.g., gfi.groovy -i ~/git/kass/input/references/genes/ -f in/MN167504_2DL1S1S2_features.fasta -g KIR2DL1S1S2 -j uniq/ -o out
 * 
 * Adding the -v option adds more output. It will include more informtion about
 * IPD names, including the resolution of the match and the version of the IPD 
 * database. GenBank will not allow this extra information.
 *
 * Input 
 *   - A fasta files with one or more sequences specific to a gene. Assumes
 *     the sequences are oriented in the 5' to 3' fashion as in IPD-KIR>
 *   - A directory of IPD-KIR reference fasta files for nuc, gen, and protein.
 *   - Gene name for the sequences (e.g., 'KIR2DL4')
 *   - A directory of files containing joint sequences that bridge features. One file for each joint.
 *   - Name of the output directory.
 *
 * Output is one table format (.tbl) as used by GenBank, etc. for each sequence.
 *
 * Requires biojava4 or later and guava.
 *   e.g., export CLASSPATH=$HOME/bin/jars/biojava4-core.jar:$CLASSPATH
 *
 * 
 * @author Dave Roe
 */
import org.biojava.nbio.core.sequence.*
import org.biojava.nbio.core.sequence.io.*
import org.biojava.nbio.core.sequence.compound.*
import org.biojava.nbio.core.sequence.transcription.*

import groovy.transform.Field
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// things that may change per run
debugging = 3 // TRACE=1, WARN=2, DEBUG=3, INFO=4, ERROR=5
@Field final Integer maxFeatureDistance = 10000 // 5000?
@Field final String NOMEN_VER = "IPD-KIR 2.9.0"
@Field final String NEW_ALLELE = "NEW"
pseudoVerbose = false // output exons for pseudogenes

// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')
File ipdDir, gffFile, outFile

OptionAccessor options = handleArgs(args)
String jointFileDir = options.j
String gene = options.g
err.println "options.v=" + options.v//todo
Boolean verbose = ((options.v != null) && (options.v != false)) ? true : false

Map<String,String> ipdProtMap // per gene
Map<String,String> ipdNucMap
Map<String,String> ipdGeneMap
(ipdProtMap, ipdNucMap, ipdGeneMap) = loadIPDMaps(new File(options.i), gene)

// gene, exon feature -> List of DNASequence
// e.g., 2DL1S1S2, exon4_3p -> [CATCACAGGTGAGAGT,...]
HashBasedTable<String,String,ArrayList<DNASequence>> jointTable = null
jointTable = loadJoints(jointFileDir, gene, fileSeparator)

// descSeqMap: description -> DNASequence unmodified from fasta file
fFile = new File(options.f)
FastaReader<DNASequence, NucleotideCompound> parentReader = new FastaReader<DNASequence, NucleotideCompound>(fFile, new PlainFastaHeaderParser<DNASequence, NucleotideCompound>(), new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet()))
LinkedHashMap<String, DNASequence> descSeqMap = parentReader.process()
if(debugging <= 3) {
    err.println "loaded ${options.f}"
    err.println "input descriptions: " + descSeqMap.keySet().join(",")
}

oFile = options.o + fileSeparator + fFile.getName().replaceFirst("_features.fasta",
                                                                 ".tbl.txt")
writer = new PrintWriter(new File(oFile).newOutputStream(), true)
descSeqMap.each { desc, dnaSeq ->
    // gene, exon feature -> zoro-based index of the start of the feature
    // e.g., 2DL1S1S2, exon4_3p -> 1227
    HashBasedTable<String,String,Integer> idxRNATable =  HashBasedTable.create()

    interpret(dnaSeq, gene, jointTable, idxRNATable)
    ArrayList<String> exonsList = idxRNATable.rowKeySet()
    ArrayList<String> sidesList = idxRNATable.columnKeySet().sort()
    if(debugging <= 3) { 
        err.println "interpreted ${exonsList.size()} exons and ${sidesList.size()} sides"
        err.println exonsList.join(',')
        err.println sidesList.join(',')
        err.println "exon0 5p = " + idxRNATable.get("exon0", "5p")
        err.println "exon1 5p = " + idxRNATable.get("exon1", "5p")
        err.println "exon10 3p = " + idxRNATable.get("exon10", "3p")
    }

    // interpret the input sequence and new annotation
    // in the context of ipd names
    ProteinSequence aaSeq
    DNASequence cdsSeq
    (idxCDSTable, aaSeq, cdsSeq, partial5p, partial3p, shortened) = joinFeatures(gene,
                                                                                 dnaSeq,
                                                                                 idxRNATable)
    // intepret utr, exon, intron
    // and label ipd names at AA and CDS
    (ipdAASet, ipdCDSSet) = interpretIPD(gene, dnaSeq, aaSeq, cdsSeq,
                                         ipdProtMap, ipdNucMap, verbose)
    // label ipd names full gene
    (fullAlleleName, bestGeneName) = interpretIPDFull(gene, dnaSeq, ipdGeneMap)
    (bestGeneName, bestAllele) = bestAllele(bestGeneName, fullAlleleName,
                                            ipdCDSSet.join("/"), ipdAASet.join("/"),
                                            verbose)
    
    output(dnaSeq, bestGeneName, bestAllele, idxRNATable, idxCDSTable, partial5p,
           partial3p, shortened, verbose, writer)
} // each sequence in input fasta
writer.close()

/* 
 * joinFeatures
 *
 * Joins the features into a full gene CDS and protein amino acids.
 *
 * @param gene String
 * @param dnaSeq DNASequence to be interpreted
 * @param idxRNATable Table of 0-based indexes for the the exons and 5p/3p
 * @return List idxCDSTable, amino acid sequence, CDS sequence, 5p partial, 3p partial, shortened string
 */
List<String> joinFeatures(String gene, DNASequence dnaSeq,
                          HashBasedTable<String,String,Integer> idxRNATable) {
    if(debugging <= 1) {
        err.println "joinFeatures(gene=${gene})"
    }
    // Initialize the Transcription Engine
    TranscriptionEngine engine = 
        new TranscriptionEngine.Builder().dnaCompounds(AmbiguityDNACompoundSet.getDNACompoundSet()).rnaCompounds(AmbiguityRNACompoundSet.getDNACompoundSet()).build()

    HashBasedTable<String,String,Integer> idxCDSTable =  HashBasedTable.create()
    
    String dnaStr = dnaSeq.getSequenceAsString()
    StringBuffer retCDS = new StringBuffer()
    String partial5p = null
    String partial3p = null
    String shortenedProtein = null
    int idxStopProtein = -1
    ProteinSequence retAA = null
    DNASequence retCDSSeq = null
    (0..10).each { exonIndex ->
        if(debugging <= 2) {
            err.println "joinFeatures: exonIndex=${exonIndex}"
        }
        String row = "exon" + exonIndex.toString()
        idx5p = idxRNATable.get(row, "5p") // e.g., exon1_5p
        idx3p = idxRNATable.get(row, "3p") // e.g., exon1_3p
        if(debugging <= 2) {
            err.println "joinFeatures: row=${row}, 5p=${idx5p}, 3p=${idx3p}"
        }

        if((idx5p == null) && (idx3p == null)) {
            if(debugging <= 2) {
                err.println "joinFeatures: initial both null"
            }

            if((partial3p == null) &&
               !((gene.contains(/3DP1/) && (exonIndex == 2)) ||
                 ((gene =~ /2DL[1-3]/) && (exonIndex == 3))  ||
                 ((gene =~ /2DS[1-5]/) && (exonIndex == 3)) ||
                 ((gene =~ /2DL2/) && (exonIndex == 3)) ||
                 (gene.contains(/2DL1S1S2/) && (exonIndex == 3)) ||
                 (gene.contains(/2DP1/) && (exonIndex == 3)) ||
                 ((gene =~ /2DL[4-5]/) && (exonIndex == 4)) ||
                 (gene.contains(/3DL3/) && (exonIndex == 6)) ||
                 (gene.contains(/3DP1/) && (exonIndex >= 6)))) {
                if(row != "exon10") {
                    // assume centromeric genes are partial on the 5'
                    // and telomeric are partial on the 3'
                    if(isTelomeric(gene)) { 
                        partial3p = previousRow(row, gene)
                    } else { 
                        partial5p = nextRow(row, gene)
                    }
                } else {
                    partial3p = "exon10"
                }
                if(debugging <= 4) {
                    err.println "joinFeatures: found partial both ${gene}, partial5p=${partial5p}, partial3p=${partial3p}"
                }
            }
            return
        } else if((idx5p == null) && (idx3p != null)) {
            if(exonIndex < 10) { 
                if(debugging <= 2) {
                    err.println "joinFeatures: found partial 5p"
                }
                idx5p = 0
                partial5p = row
            } else {
                return
            }
        } else if((idx5p != null) && (idx3p == null)) {
            if(exonIndex > 0) { 
                if(debugging <= 2) {
                    err.println "joinFeatures: found partial 3p"
                }
                idx3p = dnaStr.length() - 1
                partial3p = row
            } else {
                return
            }
        }

        if(debugging <= 2) {
            err.println "joinFeatures: idx5p=${idx5p}, idx3p=${idx3p}, partial5p=${partial5p}, partial3p=${partial3p}"
            err.println "joinFeatures: retCDS=${retCDS}"
            err.println "joinFeatures: dnaStr length=${dnaStr.length()}"
            err.println "joinFeatures: dnaStr[${idx5p}..${idx3p}]=${dnaStr[idx5p..idx3p]}"
            err.println "joinFeatures: len=${idx3p-idx5p+1}, tmpCDS=${retCDS + dnaStr[idx5p..idx3p]}"
        }

        // before adding the new full exon, translate to AA to check for
        // stop codon
        if(idxStopProtein < 0) { // if stop codon not found previously
            idxCDS5p = idx5p
            idxCDS3p = idx3p
            // test adding the full new exon
            String testStr = retCDS.toString() + dnaStr[idxCDS5p..idxCDS3p]
            DNASequence testCDSSeq = new DNASequence(testStr)
            // covert to protein
            //err.println "joinFeatures: translating ${testCDSSeq}"//todo
            retAA = engine.translate(testCDSSeq)
            retAAStr = retAA.getSequenceAsString()
            if(debugging <= 2) {
                err.println "joinFeatures: retAAStr len=${retAAStr.length()}, ${retAAStr}"
            }

            // don't check for stop codons in pseudo genes
            // or if the start of the gene is missing
            if(!(gene =~ /[23]DP/) && (partial5p == null)) {
                idxStopProtein = retAAStr.indexOf('*')
                if(debugging <= 2) {
                    err.println "joinFeatures: idxStopProtein=${idxStopProtein}"
                }

                // if you hit a stop codon before the end of the last exon
                if(idxStopProtein >= 0) {
                    idxStopProtein++  // convert to 1-based index from input
                    //old idxCDS3p = idxCDS5p + ((idxStopProtein-1) * 3)
                    idxCDS3p = idxCDS5p + ((idxStopProtein) * 3) - retCDS.length() - 1
                    if(exonIndex != 10) {
                        err.println "joinFeatures: setting shortenedProtein after stop codon, exonIndex=${exonIndex}"
                        shortenedProtein = row
                    }
                    if(debugging <= 3) {
                        err.println "joinFeatures: found stop codon, idxStopProtein=${idxStopProtein}, idxCDS5p=${idxCDS5p}, idxCDS3p=${idxCDS3p}, len=${idxCDS3p-idxCDS5p}, ${retAAStr}"
                    }
                    retAA = new ProteinSequence(retAAStr[0..idxStopProtein-1])
                    retCDSSeq = new DNASequence(dnaStr[idxCDS5p..idxCDS3p])
                }
            }
            idxCDSTable.put(row, "5p", idxCDS5p)
            idxCDSTable.put(row, "3p", idxCDS3p)
            retCDSSeq = new DNASequence(retCDS.append(dnaStr[idxCDS5p..idxCDS3p]).toString())
        }
    } // each exon

    if(debugging <= 1) {
        err.println "joinFeatures: return partial5p=${partial5p}, partial3p=${partial3p}, shortenedProtein=${shortenedProtein}"
        if(retAA != null) { 
            err.println "joinFeatures: return retAA len=${retAA.getLength()}, retCDS len=${retCDSSeq.getLength()}"
            err.println "joinFeatures: return retAA=${retAA.getSequenceAsString()}"
        }
    }
    return [idxCDSTable, retAA, retCDSSeq, partial5p, partial3p, shortenedProtein]
} // joinFeatures

// return start, end of the first 5p and last 3p of exons 0 to 10
def List<Integer> getGeneStartEnd(Table idxRNATable, Integer featureStart,
                                  Integer featureEnd) {
    if(debugging <= 1) {
        err.println "getGeneStartEnd()"
    }
    Integer retStart = null
    Integer retEnd = null

    // find start
    Integer i = new Integer(0)
    while((i <= 10) && (retStart == null)) {
        retStart = idxRNATable.get("exon" + i.toString(), "5p")
        i++
    } // start

    // find end
    i = new Integer(10)
    while((i >= 0) && (retEnd == null)) {
        retEnd = idxRNATable.get("exon" + i.toString(), "3p")
        // if 10 (UTR) is partial, set gene end to end of sequence
        if((i == 10) && (retEnd == null)) {
            retEnd = featureEnd - featureStart // added back below
            break
        }
        i--
    } // end
    retStart += featureStart + 1
    retEnd += featureStart + 1
    
    if(debugging <= 1) {
        err.println "getGeneStartEnd: return [${retStart}, ${retEnd}]"
    }
    return [retStart, retEnd]
} // getGeneStartEnd

/*
 * output
 * 
 * Write the table (tbl) file for all sequences.
 * 1 based indexes
 * http://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html
 * 
 * @param dnaSeq DNASequence to be interpreted
 * @param gene String name of the gene to which the sequence belongs
 * @param allele String name of the best allele that could be assigned
 * @param idxRNATable Table of joint mRNA sequences (exon, side) -> 0-based index
 * @param idxCDSTable Table of joint CDS sequences (exon, side) -> 0-based index
 * @param writer PrintWriter for output tbl file
 */
def void output(DNASequence dnaSeq, String gene, String allele,
                HashBasedTable<String,String,Integer> idxRNATable,
                HashBasedTable<String,String,Integer> idxCDSTable,
                String partial5p, String partial3p, String shortenedProtein,
                Boolean verbose, PrintWriter writer) {
    if(debugging <= 1) {
        err.println "output(${dnaSeq.getOriginalHeader()}, gene=${gene}, allele=${allele}, partial5p=${partial5p}, partial3p=${partial3p})"
    }

    Integer featureStart = getFeatureStart(dnaSeq.getOriginalHeader())
    Integer featureEnd = getFeatureEnd(dnaSeq.getOriginalHeader())
    if(debugging <= 2) {
        err.println "output: featureStart=${featureStart}, featureEnd=${featureEnd}"
    }
    partial5pStr = partial5p ? "<" : ""
    partial3pStr = partial3p ? ">" : ""
    partialStr = ""
    if((partial5p != null) || (partial3p != null)) {
        partialStr = " (partial)"
    }
    
    // feature header; e.g., >Feature gb|MN167504|
    writer.println ">Feature gb|${dnaSeq.getOriginalHeader()}|"

    // get start and end of gene
    Integer geneStart = null
    Integer geneEnd = null
    Integer startExonIdx = 0
    Integer endExonIdx = 10
    (geneStart, geneEnd) = getGeneStartEnd(idxRNATable, featureStart, featureEnd)

    // gene
    writer.println "${partial5pStr}${geneStart}\t${partial3pStr}${geneEnd}\tgene"
    writer.println "\t\t\tgene\t${gene}"
    writer.print "\t\t\tallele\t${gene}*${allele}"
    if(verbose == true) {
        writer.println partialStr
        writer.print "\t\t\tnote\tname checked as of ${NOMEN_VER}"
    }
    writer.println ""
    
    if((partial5p != null) || (partial3p != null)) {
        writer.println "\t\t\tnote\tpartial"
    }

    // mRNA
    if(!(gene =~ /[23]DP/)) {
        first = true
        // exon 0 is 5' utr, exon 10 is 3' utr
        idx5pStart = idxRNATable.get("exon0", "5p")
        if(idx5pStart == null) {
            idx5pStart = geneStart - featureStart - 1
        }
        idx3pEnd = idxRNATable.get("exon10", "3p")
        if(idx3pEnd == null) {
            idx3pEnd = geneEnd - featureStart - 1
        }
        (1..9).each { exonIndex ->
            String row = "exon" + exonIndex.toString()
            idx5p = idxRNATable.get(row, "5p")
            idx3p = idxRNATable.get(row, "3p")
            if(debugging <= 2) {
                err.println "output: mRNA exonIndex=${exonIndex}, row=${row}, idx5p=${idx5p}, idx3p=${idx3p}"
            }
            if((idx5p == null) && (idx3p == null)) {
                return // next exon
            }
            idx5pNew = idx5p ? (idx5p + featureStart + 1) : null
            idx3pNew = idx3p ? (idx3p + featureStart + 1) : null
            String idx5pStr = idx5p ? idx5pNew : ""
            String idx3pStr = idx3p ? idx3pNew : ""
            if(row == "exon1") { // extend 5' through utr
                idx5pStr = idx5pStart + featureStart + 1
            } else if(row == "exon9") { // extend 3' through utr
                idx3pStr = idx3pEnd + featureStart + 1
            }
            partial5pStr = (partial5p == row) ? '<' : ''
            partial3pStr = (partial3p == row) ? '>' : ''
            if((row == "exon9") && (partial3p == "exon10")) {
                partial3pStr = '>'
            }

            if(debugging <= 3) {
                err.println "output: (partial5p == row)"
                err.println "output: partial5p=${partial5p}, partial5pStr=${partial5pStr}"
                err.println "output: partial3p=${partial3p}, partial3pStr=${partial3pStr}"
            }
            
            sectionStr = ""
            if(first == true) { 
                sectionStr = "\tmRNA"
                first = false
            }
            writer.println  "${partial5pStr}${idx5pStr}\t${partial3pStr}${idx3pStr}" + sectionStr
        } // each mRNA exon
        writer.println "\t\t\tproduct\tkiller cell immunoglobulin-like receptor"
    } // pseudo

    // CDS
    first = true
    pseudoExonCount = 1
    (1..9).each { exonIndex ->
        if(debugging <= 2) {
            err.println "output: CDS exonIndex=${exonIndex}"
        }
        String row = "exon" + exonIndex.toString()
        idx5p = idxCDSTable.get(row, "5p")
        idx3p = idxCDSTable.get(row, "3p")
        if((idx5p == null) && (idx3p == null)) {
            return // next exon
        }
        idx5pNew = idx5p ? (idx5p + featureStart + 1) : null
        idx3pNew = idx3p ? (idx3p + featureStart + 1) : null
        idx5pStr = idx5p ? idx5pNew : ""
        idx3pStr = idx3p ? idx3pNew : ""
        partial5pStr = (partial5p == row) ? '<' : ''
        partial3pStr = (partial3p == row) ? '>' : ''

        sectionStr = ""
        if(first == true) { 
            sectionStr = "\tCDS"
            first = false
        }
        // exon for pseudo genes
        if(gene =~ /[23]DP1/) {
            sectionStr = "\texon"
        }

        if(!(gene =~ /[23]DP1/) ||
           ((gene =~ /[23]DP1/) && (pseudoVerbose == true))) {
            writer.println  "${partial5pStr}${idx5pStr}\t${partial3pStr}${idx3pStr}" +
                sectionStr
        }
        // exon number for pseudo genes
        if((gene =~ /[23]DP1/) && (pseudoVerbose == true)) {
            writer.println "\t\t\tnumber\t${pseudoExonCount}"
            pseudoExonCount++
        }
    } // each CDS exon
    if(!(gene =~ /[23]DP1/)) {
        writer.println "\t\t\tproduct\tkiller cell immunoglobulin-like receptor"
        writer.println "\t\t\ttransl_table\t1"
    } else {
        writer.println "\t\t\tgene_desc\tkiller cell immunoglobulin-like receptor"
        writer.println "\t\t\tpseudo"
        writer.println "\t\t\tpseudogene\tunprocessed"
    }
    if(debugging <= 1) {
        err.println "output: return"
    }
    return
} // output

/*
 * interpret
 * Documents the locations of gene features in the given sequence.
 * 
 * @param dnaSeq DNASequence to be interpreted
 * @param gene String name of the gene to which the sequence belongs
 * @param jointTable Table of joint sequences (exon, side) -> List<DNASequences>
 * @param idxRNATable Table of joint sequences (exon, side) -> 0-based index
 */
def void interpret(DNASequence dnaSeq, String gene, 
                   HashBasedTable<String,String,ArrayList<DNASequence>> jointTable,
                   HashBasedTable<String,String,Integer> idxRNATable) {
    if(debugging <= 1) {
        err.println "interpret(${dnaSeq.getOriginalHeader()}, ${gene})"
        err.println "interpret: dnaSeq length=${dnaSeq.getLength()}"
    }

    dnaStr = dnaSeq.getSequenceAsString()
    //err.println dnaStr //todo
    Integer lastIdx = 0
    (0..10).each { exonIndex ->
        if(debugging <= 2) {
            err.println "interpret: exonIndex=${exonIndex}"
        }
        ["5p", "3p"].each { side ->
            String row = "exon" + exonIndex.toString()
            //err.println "row=${row}|, col=${side}|"//todo
            //err.println jointTable.get("exon0", "5p")//todo
            queryList = jointTable.get(row, side) // e.g., exon1_5p
            if(queryList == null) {
                if(debugging <= 3) {
                    err.println "interpret: no index found for ${row} ${side}"
                }
                return
            }
            if(debugging <= 3) {
                err.println "interpret: queryList=" + queryList.join(',')
            }

            // only search toward the 3p direction
            // shorten the string with every find
            Integer idx = find(dnaStr[lastIdx..-1], queryList, side)
            if(idx >= 0) { 
                if(debugging <= 3) {
                    err.println "interpret: found exon${exonIndex}, ${side} for ${dnaSeq.getOriginalHeader()}: idx(len-1)=$idx"
                }
                Integer setVal = new Integer(lastIdx + idx)
                idxRNATable.put(row, side, setVal)
                if(debugging <= 2) {
                    err.println "row=${row}, setVal=${setVal}"
                    err.println "checking " + "exon${exonIndex}:" + idxRNATable.get("exon${exonIndex}", side)
                    err.println "checking exon0: " + idxRNATable.get("exon0", side)
                }

                lastIdx += idx
            } else {
                // this isn't always an error; e.g., 3DP1 exon 2
                if(debugging <= 5) {
                    err.println "WARNING: couldn't find ${row}, ${side} for ${dnaSeq.getOriginalHeader()}"
                    err.println queryList.join(',')
                }
            }
        }
    } // each exon and utr

    if(debugging <= 1) {
        err.println "interpret: return"
    }
} // interpret

/*
 * find
 * Finds one of the query sequences in the target and returns the 0-based index
 * _of the junction point_.
 * Assumes features cannot be more than maxFeatureDistance bp away.
 *
 * @param String target sequence
 * @param String query sequence
 * @param String side 5p or 3p
 * @return int of the 0-based index, or < 0 if none found
 */
int find(String target, ArrayList<String> queryList, String side) {
    if(debugging <= 1) {
        err.println "find(target=${target}, side=${side})"
    }
    int ret = -1
    targetLen = target.length()

    for(Iterator qiter = queryList.iterator(); qiter.hasNext();) {
        DNASequence query = qiter.next()
        queryStr = query.getSequenceAsString()
        queryLen = queryStr.length()
        // assume feature cannot be more than maxFeatureDistance bp away
        end = maxFeatureDistance
        if(targetLen - 1 < maxFeatureDistance) {
            end = targetLen - 1
        }
        ret = target[0..end].toUpperCase().indexOf(queryStr.toUpperCase())
        if(debugging <= 3) {
            err.println "targetLen=${targetLen}, end=${end}"
            err.println "find: ${queryStr} returned ${ret}"
        }
        if(ret >= 0) {
            // count from intron side when switching to inexact matches (todo)
            ret += queryLen/2 // 5p: start of the feature
            if(side == "3p") { // end of the feature
                ret--
            }
            if(debugging <= 3) {
                err.println "find: set ${queryStr} to ${ret}"
            }
            break
        }
    } // each query sequence

    if(debugging <= 1) {
        err.println "find: return $ret"
    }
    return ret
} // find

/* 
 * Checks a sequence to see if it matches an ipd-kir 'nuc' sequence.
 *
 * @return List with 3 values:  ipdGeneName(String), ccall(String), fullName(String)
 */
def ArrayList<String> getcDNAMatches(DNASequence cSeq, Map<String, String> ipdNucMap) { 
    if(debugging <= 1) {
        err.println "getcDNAMatches(cSeq=${cSeq})"
    }
    String ipdGeneNew
    String ccall
    String fullName
    ipdNucMap.each { nucHeader, nucSeq ->
        /* Use starts-with instead of equals because the 
         IPD cDNA sequences are some times too long.
         e.g., KIR2DL4*0080203
         */
        if(sequenceContains(cSeq, nucSeq)) { // cDNA match *
            if(debugging <= 2) {
                err.println "getcDNAMatches: found *cDNA match: ${nucHeader}"
            }
            (ipdGeneNew, ccall) = truncateToCDNA(nucHeader)
            fullName = "${ipdGeneNew}*${ccall}"
        } // if cDNA match
    } // each IPD/IMGA _nuc sequence
    if(debugging <= 1) {
        err.println "getcDNAMatches: return: [ipdGeneNew=${ipdGeneNew}, ccall=${ccall}, fullName=${fullName}]"
    }
    return [ipdGeneNew, ccall, fullName]
} // getcDNAMatches

/*
 * Parses and interprets a sequence wrt IPD-KIR names.
 *
 * @param gene String gene name
 * @param ipdProtMap the IPD _prot.fasta sequences
 * @param ipdNucMap the IPD _nuc.fasta sequences
 * @param ipdGeneMap the IPD _gen.fasta sequences
 * @return List<String> aaSet and cdsSet containing ipd names that match
 */
def List interpretIPD(String gene, DNASequence dnaSeq, ProteinSequence aaSeq,
                      DNASequence cdsSeq, Map<String, String> ipdProtMap,
                      Map<String, String> ipdNucMap, Boolean verbose) {
    if(debugging <= 1) {
        err.println "interpretIPD(gene=${gene}, dnaSeq=${dnaSeq.getOriginalHeader()})"
    }
    desc = dnaSeq.getOriginalHeader()
    // return these
    HashSet<String> pCallSet = new HashSet()
    HashSet<String> cCallSet = new HashSet()
    ipdProtMap.each { header, seq ->
        if(debugging <= 1) {
            err.println "interpretIPD: protein ${header}=${seq}"
        }
        protTrue = sequenceEquals(aaSeq, seq)
        ipdGene = gene
        // protein match *
        if(protTrue) {
            // check cDNA sequence
            ipdGeneNew = ipdGene
            if(debugging <= 2) { 
                err.println "interpretIPD: ipdGene 1 =${ipdGene}"
            }
            (ipdGeneNew, pcall) = truncateToProtein(header)
            gene = ipdGeneNew
            fullName = "${ipdGeneNew}*${pcall}"
            if(debugging <= 3) {
                err.println "interpretIPD: *protein match: ${fullName}"
            }
            if(debugging <= 1) {
                err.println "interpretIPD: dnaSeq=${dnaSeq}"
                err.println "interpretIPD: ipd seq=${seq}"
            }
            pCallSet.add(fullName)
            // find the cDNA sequence
            (ipdGeneRet, ccall, fullName) = getcDNAMatches(cdsSeq, ipdNucMap)
            if(ipdGeneRet != null) {
                gene = ipdGeneRet
                cCallSet.add(fullName)
                if(debugging <= 1) { 
                    err.println "interpretIPD: adding ${fullName} to ${desc} ccall set, size now ${cCallSet.size()}"
                }
            }
        } // protein match
    } // each ipd protein sequence for this gene

    if(pCallSet.size() == 0) {
        pCallSet.add("${gene}*" + NEW_ALLELE)
        
    }
    if(cCallSet.size() == 0) {
        cCallSet.add("${gene}*" + NEW_ALLELE)
    }
    if(debugging <= 1) {
        err.println "interpretIPD: return pCallSet=${pCallSet.join(',')}, cCallSet=${cCallSet.join(',')}"
    }
    return [pCallSet, cCallSet]
} // interpretIPD

/*
 * Annotates the 'full gene' sequence wrt IPR alleles.
 *
 * @param gene String gene name
 * @param dnaSeq DNASequence the sequence being interpreted
 * @param ipdGeneMap Map of the allele names to full gene alleles
 * @param cCallSet the cDNA call(s) for this gene
 * @param pCallSet the protein call(s) for this gene
 * @return Set of allele names that match at full allele (or partial full allele) and new gene name
 */
def List interpretIPDFull(String gene, DNASequence dnaSeq,
                          Map<String,String> ipdGeneMap) {
    if(debugging <= 1) {
        err.println "interpretIPDFull(gene=${gene}, dnaSeq=${dnaSeq.getOriginalHeader()})"
    }
    genSeqStr = dnaSeq.getSequenceAsString()
    if(debugging <= 3) {        
        err.println "interpretIPDFull: ${gene} sequence is " +
            "${dnaSeq.getLength()} bp long"
    }
    
    ipdGene = gene
    HashSet<String> gCallSet = new HashSet()
    String gcall = null
    ipdGeneMap.each { header, seq ->
        // full sequence match *
        if(sequenceContains(seq, dnaSeq)) {
            (ipdGeneNew, gcall) = interpFullGene(header, gene)
            fullName = "${ipdGeneNew}*${gcall}"
            if(debugging <= 3) {
                err.println "interpretIPDFull: *match on full gene: ${fullName}"
            }
            gCallSet.add(fullName)
            (cGene, ccall) = truncateToCDNA(fullName)
            fullCDNA = "${cGene}*${ccall}"
            //cCallSet.add(fullCDNA)
            (ipdGeneNew, pcall) = truncateToProtein(fullCDNA)
            fullProt = "${ipdGeneNew}*${pcall}"
            //pCallSet.add(fullProt)
            ipdGene = ipdGeneNew
        }
    } // each ipd/imgt full sequence

    if(!gCallSet.isEmpty()) {
        gcall = gCallSet.sort().join("/")
    } else {
        gcall = "${ipdGene}*" + NEW_ALLELE
    }
    
    if(debugging <= 1) {
        err.println "interpretIPDFull: return ${gcall}"
    }
    return [gcall, ipdGene]
} // interpretIPDFull

/* 
 * loadJoints
 * Fills a Map with sequences representing the bridge or joint between 
 * two features.
 *
 * @param dir String of directory to the files containing the sequences
 * @param gene String name of the gene
 * @return Table of the joint sequences;  (e.g., exon1, 3p) -> List of DNASequences
 */
HashBasedTable<String,String,ArrayList<DNASequence>> loadJoints(String dir,
                                                                String gene,
                                                                String fileSeparator) {
    if(debugging <= 1) {
        err.println "loadJoints(dir=${dir}, gene=${gene})"
    }
    HashBasedTable<String,String,ArrayList<DNASequence>> retTable = HashBasedTable.create()

    jGene = gene.replaceFirst("KIR", "")
    // 2DL2L3S3S4S5
    jGene = jGene.replaceFirst("2DS4", "2DL2L3S3S4S5")
    
    //pName = dir + fileSeparator + jGene + "*\\.txt"
    // todo: make non-recursive
    lFiles = new FileNameByRegexFinder().getFileNames(dir, /${jGene}/)
    if(debugging <= 2) {
        err.println "jGene=${jGene}"
        err.println "loadJoints: found ${lFiles}"
    }
    if(lFiles.size() == 0) {
        err.println "ERROR: couldn't find ${jGene} joint files in ${dir}"
        System.exit(1)
    }
    first = true
    lFiles.each { f ->
        // name without the extension (e.g., 2DL1S1S2_exon4_3p.txt)
        eName = f.replaceFirst(~/\.[^\.]+$/, '')
        (geneStr, exonStr, sideStr) = eName.split('_')
        seqList = new ArrayList()

        new File(f).eachLine { line ->
            lineT = line.trim()
            if((lineT.length() % 2) != 0) {
                err.println "ERROR: joint sequence has uneven length: ${lineT}, ${f}"
                System.exit(1)
            }
            seqList.add(new DNASequence(lineT))
        } // each line
        if(debugging <= 3) {
            err.println "loadJoints: adding ${exonStr}, ${sideStr} = ${seqList.size()} sequences"
        }
        retTable.put(exonStr, sideStr, seqList)
    } // each file

    if(debugging <= 1) {
        err.println "loadJoints: return"
    }
    return retTable
} // loadJoints

/*
 * Creates three maps for the '_prot', '_nuc', and '_gen' IPD files.
 * Each map is allele name -> sequence.
 *
 * @param ipdDir File of the directory of the '_prot.fasta', '_nuc.fasta',
 *               and '_gen.fasta' files
 * @param ipdGene String of the IPD gene name
 * @return ArrayList<Map> three for prot, nuc, and gen.
 *                        prot: LinkedHashMap<String, ProteinSequence>
 *                        nuc & gen: LinkedHashMap<String, DNASequence>
 *                        Each is fasta header -> sequence
 */
def ArrayList<Map>loadIPDMaps(File ipdDir, String ipdGene) {
    if(debugging <= 1) {
        err.println "loadIPDMaps(ipdDir=${ipdDir}, ipdGene=${ipdGene})"
    }
    // the protein (amino acids) sequences (_prot.fasta)
    protFileName = ipdDir.getAbsolutePath() + "/${ipdGene}_prot.fasta"
    LinkedHashMap<String, ProteinSequence> protSeqMap =
        FastaReaderHelper.readFastaProteinSequence(new File(protFileName))
    if(debugging <= 3) {        
        err.println protSeqMap.keySet().size() + " proteins in ${protFileName}"
    }

    // the nucleic acid (dna) sequences for the coding regions/exons (_nuc.fasta)
    nucFileName = ipdDir.getAbsolutePath() + "/${ipdGene}_nuc.fasta"
    LinkedHashMap<String, DNASequence> nucSeqMap =
        FastaReaderHelper.readFastaDNASequence(new File(nucFileName))
    if(debugging <= 3) {        
        err.println nucSeqMap.keySet().size() + " nuc sequences in ${nucFileName}"
    }

    // the full nucleic acid (dna) sequences (_gen.fasta)
    genFileName = ipdDir.getAbsolutePath() + "/${ipdGene}_gen.fasta"
    LinkedHashMap<String, DNASequence> genSeqMap =
        FastaReaderHelper.readFastaDNASequence(new File(genFileName))
    if(debugging <= 3) {        
        err.println genSeqMap.keySet().size() + " gene sequences in ${genFileName}"
    }

    if(debugging <= 1) {
        err.println "loadIPDMaps: return"
    }
    return [protSeqMap, nucSeqMap, genSeqMap]
} // loadIPDMaps

// todo: document
def Boolean sequenceEquals(org.biojava.nbio.core.sequence.template.Sequence a,
                           org.biojava.nbio.core.sequence.template.Sequence b) {
    if(debugging <= 1) {
        err.println "sequenceEquals()"
    }
    Boolean ret = false
    if((a == null) && (b == null)) {
        return true
    } else if(((a != null) && (b == null)) || ((a == null) && (b != null))) {
        return false
    }
    // augustus puts a '*' at the end of proteins
    astr = a.getSequenceAsString().replaceAll("\\*", "")
    bstr = b.getSequenceAsString().replaceAll("\\*", "")
    bstrrc = bstr
    if(b.class =~ /Protein/) { // protein, not DNA (or RNA)
        // IPD-KIR sequences have trailing 'X's; remove them if present
        astr = astr - ~/X$/
        bstr = bstr - ~/X$/
    } else { // dna
        bstrrc = b.getReverseComplement().getSequenceAsString()
    }

    if(astr.equalsIgnoreCase(bstr) || astr.equalsIgnoreCase(bstrrc)) {
        if(debugging <= 1) {
            err.println "sequenceEquals: found match"
        }
        ret = true
    }
    if(debugging <= 1) {
        err.println "sequenceEquals: return ${ret}"
    }
    return ret
} // sequenceEquals

// todo: document
def Boolean sequenceStartsWith(org.biojava.nbio.core.sequence.template.Sequence a,
                               org.biojava.nbio.core.sequence.template.Sequence b) {
    Boolean ret = false
    astr = a.getSequenceAsString()
    bstr = b.getSequenceAsString()
    bstrrc = bstr
    if(b.class =~ /Protein/) { // protein, not DNA (or RNA)
        // IPD-KIR sequences have trailing 'X's; remove them if present
        astr = astr - ~/X$/
        bstr = bstr - ~/X$/
    } else { // dna
        bstrrc = b.getReverseComplement().getSequenceAsString()
    }

    if(astr.toLowerCase().startsWith(bstr.toLowerCase()) ||
       astr.toLowerCase().startsWith(bstrrc.toLowerCase()) ||
       bstr.toLowerCase().startsWith(astr.toLowerCase()) ||
       bstrrc.toLowerCase().startsWith(astr.toLowerCase())) {
        ret = true
    }
    return ret
} // sequenceStartsWith

// todo: document
def Boolean sequenceContains(org.biojava.nbio.core.sequence.template.Sequence a,
                             org.biojava.nbio.core.sequence.template.Sequence b) {
    Boolean ret = false
    astr = a.getSequenceAsString()
    bstr = b.getSequenceAsString()
    bstrrc = bstr
    if(b.class =~ /Protein/) { // protein, not DNA (or RNA)
        // IPD-KIR sequences have trailing 'X's; remove them if present
        astr = astr - ~/X$/
        bstr = bstr - ~/X$/
    } else { // dna
        bstrrc = b.getReverseComplement().getSequenceAsString()
    }

    if(astr.toLowerCase().contains(bstr.toLowerCase()) ||
       astr.toLowerCase().contains(bstrrc.toLowerCase()) ||
       bstr.toLowerCase().contains(astr.toLowerCase()) ||
       bstrrc.toLowerCase().contains(astr.toLowerCase())) {
        if(debugging <= 1) { 
            err.println "sequenceContains: true for ${astr} and ${bstr}"
        }
        ret = true
    }
    return ret
} // sequenceContains

/*
 * Returns the protein version of the given KIR or HLA allele.
 *
 * @param alelle a String containing a '*'-separated allele name 
 *               (e.g., KIR2DL4*0010201, C*01:02:05)
 * @return a List containing two Strings: the locus name from the fasta file and
 *         a String containing just the post-separator allele contents
 *         at the cDNA level (e.g., 0010201, 01:02:05). 
 *         The gene name may be different than the input
 *         name: e.g., 2DL5 in, 2DL5A or 2DL5B out. Returns null Strings
 *         if not found.
 */
def ArrayList<String> truncateToCDNA(String allele) {
    if(debugging <= 1) {
        err.println "truncateToCDNA(allele=${allele})"
    }
    //patternStr = "(.*)\\*(.*)"
    patternStr = "(.*)_(.*)"
    if(allele =~ /KIR2DL5/) {
      // problem here with the "KIR" part
      patternStr = patternStr.replace("2DL5", "2DL5[A-B]")
    }
    pattern = ~/${patternStr}/
    if(debugging <= 2) {
        err.println "truncateToCDNA: ipd allele matching " +
            "${pattern} against ${allele}"
    }
    matcher = pattern.matcher(allele)
    if(matcher.size() == 0) {
        // e.g. KIR3DL2*018
        patternStr = "(.*)\\*(.*)"
        if(allele =~ /KIR2DL5/) {
            patternStr = patternStr.replace("KIR2DL5", "KIR2DL5[A-B]")
        }
        pattern = ~/${patternStr}/
        if(debugging <= 2) {
            err.println "truncateToProtein: ipd allele matching " +
                "${pattern} against ${allele}"
        }
        matcher = pattern.matcher(allele)
        if(matcher.size() == 0) {

            if(debugging <= 3) {
                err.println "truncateToCDNA: ipd allele matching " +
                    "${pattern} against ${allele}"
                err.println "truncateToCDNA: return null, nothing found"
            }
            return [null, null]
        }
    }
    locus = matcher[0][1]
    name = matcher[0][2]

    if(debugging <= 2) {
        err.println "truncateToCDNA: locus=${locus}, name=${name}"
    }
 
    if(allele =~ /^KIR[2-3]/) {
        // e.g., 0090102 -> 00901
        nl = name.length()
        if(nl > 5) {
            nl = 5
        }
        allele = name.substring(0, nl)
     } else {
        // e.g., 17:01:01:03 -> 17:01:01
        fields = name.split(':')
        allele = "${fields[0]}:${fields[1]}"
        if(fields[2] != null) {
            allele += ":${fields[2]}"
        }
    }
    if(debugging <= 1) {
        err.println "truncateToCDNA: return ${allele}"
    }
    return [locus, allele]
} // truncateToCDNA

/*
 * Returns the protein version of the given KIR or HLA allele.
 *
 * http://www.ebi.ac.uk/ipd/kir/alleles.html
 * http://hla.alleles.org/nomenclature/naming.html
 * 
 * @param header a String containing a '*'-separated allele name 
 *               (e.g., KIR2DL4*0010201, C*01:02:05); can contain other stuff
 * @return a List of Strings containing the IPD gene name and just the 
 *         post-separator allele contents at the protein level 
 *         (e.g., 001, 01:02). The gene name may be different than the input
 *         name: e.g., 2DL5 in, 2DL5A or 2DL5B out. Returns null Strings
 *         if not found.
 */
def ArrayList<String> truncateToProtein(String allele) {
    if(debugging <= 1) {
        err.println "truncateToProtein(allele=${allele})"
    }
    // e.g. KIR2DP1_00101
    patternStr = "(.*)_(.*)"

    if(allele =~ /KIR2DL5/) {
        patternStr = patternStr.replace("KIR2DL5", "KIR2DL5[A-B]")
    }
    pattern = ~/${patternStr}/
    if(debugging <= 2) {
        err.println "truncateToProtein: ipd allele matching " +
            "${pattern} against ${allele}"
    }
    matcher = pattern.matcher(allele)
    if(matcher.size() == 0) {
        // e.g. KIR2DP1*00101
        patternStr = "(.*)\\*(.*)"
        if(allele =~ /KIR2DL5/) {
            patternStr = patternStr.replace("KIR2DL5", "KIR2DL5[A-B]")
        }
        pattern = ~/${patternStr}/
        if(debugging <= 2) {
            err.println "truncateToProtein: ipd allele matching " +
                "${pattern} against ${allele}"
        }
        matcher = pattern.matcher(allele)
        if(matcher.size() == 0) {
            if(debugging <= 3) {
                err.println "truncateToProtein: ipd allele matching " +
                    "${pattern} against ${allele}"
                err.println "truncateToProtein: return null, nothing found"
            }
            return [null, null]
        }
    }
//    accession = matcher[0][1]
    locus = matcher[0][1]
    name = matcher[0][2]

    if(debugging <= 2) {
        err.println "truncateToProtein: locus=${locus}, name=${name}"
    }
 
    if(allele =~ /^KIR[2-3]/) { // KIR check
        // e.g., 0090102 -> 3DP1*009
        allele = name.substring(0, 3)
     } else {
        // e.g., 78:01:01 -> B*78:01
        fields = name.split(':')
        allele = "${fields[0]}:${fields[1]}"
    }
    if(debugging <= 1) {
        err.println "truncateToProtein: return [$locus, ]${allele}]"
    }
    return [locus, allele]
} // truncateToProtein

/*
 * Returns the full allele name given fasta header.
 * 
 * For KIR, the 'KIR' is not in the allele names in the *_gen.fasta files.
 * e.g., 'KIR2DL4' in _nuc.fasta and _prot.fasta is '2DL4' in _gen.fasta.
 *
 * @param allele locus_allele (eg., KIR2DS2_0010101)
 * @return a List containing two Strings: the locus name from the fasta file and
 *         a String containing the full allele name if found. The gene name may be 
 *         different than the input name: e.g., 2DL5 in, 2DL5A or 2DL5B out. 
 *         Returns null Strings if not found.
 */
def ArrayList<String> interpFullGene(String allele, String gene) {
    if(debugging <= 1) {
        err.println "interpFullGene(allele=${allele})"
    }

    if(gene =~ /KIR/) {
        gene = gene.replace("KIR", "")
        if(gene == "KIR2DL5") {
            gene = gene.replace("2DL5", "2DL5[A-B]")
        }
    }
    patternStr = "(.*)_(.*)"
    pattern = ~/${patternStr}/
    if(debugging <= 2) {
        err.println "interpFullGene: ipd allele matching " +
            "${pattern} against ${allele}"
    }
    matcher = pattern.matcher(allele)
    if(matcher.size() == 0) {
        if(debugging <= 3) {
            err.println "interpFullGene: ipd allele matching " +
                "${pattern} against ${allele}"
            err.println "interpFullGene: return null, nothing found"
        }
        return [null, null]
    }
//    accession = matcher[0][1]
    locus = matcher[0][1]
    name = matcher[0][2]

    if(debugging <= 2) {
        err.println "interpFullGene: locus=${locus}, name=${name}"
    }
    /*if(!locus.contains("KIR")) {
        locus = "KIR${locus}"
    }*/
    if(debugging <= 1) {
        err.println "interpFullGene: return ${allele}"
    }
    return [locus, name]
} // interpFullGene

/*
 * bestAllele
 * 
 * @return List<String> best gene name and best allele name
 */
def List<String> bestAllele(String gene, String gcall, String ccall, String pcall,
                            Boolean verbose) {
    // get best IPD-KIR call
    String bestGene = null
    String bestAllele = null
    String bestRes = null
    gcall += (verbose == true) ? " (full gene)" : ""
    ccall += (verbose == true) ? " (cDNA)" : ""
    pcall += (verbose == true) ? " (protein)" : ""
    if(!gcall.contains("NEW")) {
        bestAllele = gcall
        bestRes = "full gene"
    } else if(!ccall.contains("NEW")) {
        bestAllele = ccall
        bestRes = "cDNA"
//    } else if(!pcall.contains("NEW")) {
    } else {      
        bestAllele = pcall
        bestRes = "protein"
    }


    if(debugging <= 1) {
        err.println "interpFullGene: bestAllele=${bestAllele}"
    }    
    if(bestAllele != null) { 
        (bestGene, allele) = bestAllele.split('\\*')
    }
    return [bestGene, allele]
} // bestAllele

def Boolean isTelomeric(String gene) {
    tList = ["3DL2", "2DS1", "2DS3", "2DS4", "2DS5", "2DL5", "3DL1", "3DS1", "2DL4"]
    ret = new Boolean(false)
    for(Iterator giter = tList.iterator(); giter.hasNext();) {
        String g = (String)giter.next()
        if(gene.contains(g)) {
            ret = new Boolean(true)
            break
        }
    }
    return ret
} // isTelomeric

def String previousRow(String row, String gene) {
    if(debugging <= 1) {
        err.println "previousRow(row=${row}, gene=${gene})"
    }
    Integer exonIndex = new Integer(row[-1])
    lastTwoChars = row[-2..-1] // check for "exon10"
    if(lastTwoChars == "10") {
        exonIndex = new Integer(lastTwoChars)
    }
    exonIndex--
    if((gene.contains(/3DP1/) && (exonIndex == 2)) ||
       ((gene =~ /2DL[1-3]/) && (exonIndex == 3))  ||
       ((gene =~ /2DS[1-5]/) && (exonIndex == 3)) ||
       ((gene =~ /2DL2/) && (exonIndex == 3)) ||
       (gene.contains(/2DL1S1S2/) && (exonIndex == 3)) ||
       (gene.contains(/2DP1/) && (exonIndex == 3)) ||
       ((gene =~ /2DL[4-5]/) && (exonIndex == 4)) ||
       (gene.contains(/3DL3/) && (exonIndex == 6)) ||
       (gene.contains(/3DP1/) && (exonIndex >= 6))) { 
        exonIndex--
    }
    ret = "exon" + exonIndex.toString()
    if(debugging <= 1) {
        err.println "previousRow: ret=${ret}"
    }
    return ret
} // previousRow

def String nextRow(String row, String gene) {
    if(debugging <= 1) {
        err.println "nextRow(row=${row}, gene=${gene})"
    }

    Integer exonIndex = new Integer(row[-1])
    lastTwoChars = row[-2..-1] // check for "exon10"
    if(lastTwoChars == "10") {
        return null
    }
    exonIndex++
    if((gene.contains(/3DP1/) && (exonIndex == 2)) ||
       ((gene =~ /2DL[1-3]/) && (exonIndex == 3))  ||
       ((gene =~ /2DS[1-5]/) && (exonIndex == 3)) ||
       ((gene =~ /2DL2/) && (exonIndex == 3)) ||
       (gene.contains(/2DL1S1S2/) && (exonIndex == 3)) ||
       (gene.contains(/2DP1/) && (exonIndex == 3)) ||
       ((gene =~ /2DL[4-5]/) && (exonIndex == 4)) ||
       (gene.contains(/3DL3/) && (exonIndex == 6)) ||
       (gene.contains(/3DP1/) && (exonIndex >= 6))) { 
        exonIndex++
    }
    ret = "exon" + exonIndex.toString()
    if(debugging <= 1) {
        err.println "nextRow: ret=${ret}"
    }
    return ret
} // nextRow

// MN167504_2DL2L3S3S4S5_23918-43868
def Integer getFeatureStart(String desc) {
    if(debugging <= 3) {
        err.println "getFeatureStart(desc=${desc})"
    }

    // separate location
    idx = desc.lastIndexOf('_')
    posTmp = desc[idx+1..-1]
    (posTmp, end) = posTmp.split('-')
    Integer position = posTmp.toInteger()
    if(debugging <= 1) {
        err.println "getFeatureStart: return ${position}"
    }
    return position
} // getFeatureStart

// MN167504_2DL2L3S3S4S5_23918-43868
def Integer getFeatureEnd(String desc) {
    if(debugging <= 3) {
        err.println "getFeatureEnd(desc=${desc})"
    }

    // separate location
    idx = desc.lastIndexOf('_')
    posTmp = desc[idx+1..-1]
    (start, posTmp) = posTmp.split('-')
    Integer position = posTmp.toInteger()
    if(debugging <= 1) {
        err.println "getFeatureEnd: return ${position}"
    }
    return position
} // getFeatureEnd

/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'gfi.groovy [options] ',
      header:'Options:')
    cli.help('print this message')
    cli.i(longOpt:'ipd', args:1, argName:'ipd', 'input IPD directory',
      required: true)
    cli.f(longOpt:'fasta', args:1, argName:'fasta', 'input fasta file',
      required: true)
    cli.g(longOpt:'gene', args:1, argName:'gene', "gene name (e.g., 'KIR2DL4')",
      required: true)
    cli.j(longOpt:'joints', args:1, argName:'joints', "input directory for files of joint sequences",
      required: true)
    cli.v(longOpt:'verbose', args:1, argName:'verbose', "verbose output",
      required: false)
    cli.o(longOpt:'out', args:1, argName:'output', 'output directory',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
