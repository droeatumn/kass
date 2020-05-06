#!/usr/bin/env groovy
/*
 * Convert AUGUSTUS GFF output to NCBI's table2asn_GFF. Also interpret the sequence
 * in the context of IPD-KIR allele names.
 *
 * Requires 
 *   - Commandline: http://www.dishevelled.org/commandline/
 *   - Super CSV: http://super-csv.github.io/super-csv/
 *   - BioJava 4's Core jar: http://www.biojava/org
 *   - slf4j-api & slf4j-log4j12 http://www.slf4j.org/
 *   - log4j 1 https://logging.apache.org/log4j/1.2/
 *   e.g., 
 *     export CLASSPATH=$HOME/bin/jars/dsh-commandline-1.1.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/biojava4-core.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/dsh-commandline-1.1.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/super-csv.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/biojava4-core.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/slf4j-api-1.7.5.jar:$CLASSPATH
 * @author Dave Roe
 */

import org.biojava.nbio.core.sequence.*
import org.biojava.nbio.core.sequence.io.*
import org.biojava.nbio.core.sequence.compound.*

import org.supercsv.io.*
import org.supercsv.prefs.*
import org.supercsv.cellprocessor.ift.*

import org.dishevelled.commandline.*
import org.dishevelled.commandline.argument.*

import groovy.transform.Field

// things that may change per run
debugging = 3 // TRACE=1, WARN=2, DEBUG=3, INFO=4, ERROR=5

@Field final String GENE_SYSTEM
@Field final String GENE_NAME
@Field final String KIR_NOMEN_VER = "IPD-KIR 3.9.0"
@Field final String HLA_NOMEN_VER = "IMGT/HLA 3.22.0"
@Field String NOMEN_VER // tbd by input: KIR_NOMEN_VER or HLA_NOMEN_VER
// for now, just all the possible sizes
@Field final KIR_EXON_SIZES = [34, 40, 36, 285, 300, 299, 294, 51, 102, 105, 53,
                               177, 42, 270, 8, 126, 210]
@Field final String NEW_ALLELE = "NEW"

err = System.err
File ipdDir, gffFile, outFile

// startIndex is a 1-based index
(ipdDir, gffFile, seqFile, outGLFile, GENE_NAME) = handleArgs(args)

if(GENE_NAME.contains("KIR")) {
    NOMEN_VER = KIR_NOMEN_VER
    GENE_SYSTEM = "KIR"
} else {
    NOMEN_VER = HLA_NOMEN_VER
    GENE_SYSTEM = "HLA"
}

Map<String,String> ipdProtMap // per gene
Map<String,String> ipdNucMap
Map<String,String> ipdGeneMap
(ipdProtMap, ipdNucMap, ipdGeneMap) = loadIPDMaps(ipdDir, GENE_NAME)

// description -> DNASequence unmodified from fasta file
FastaReader<DNASequence, NucleotideCompound> parentReader = new FastaReader<DNASequence, NucleotideCompound>(seqFile, new PlainFastaHeaderParser<DNASequence, NucleotideCompound>(), new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet()))
LinkedHashMap<String, DNASequence> descSeqMap = parentReader.process()
if(debugging <= 3) { 
    err.println "input descriptions: " + descSeqMap.keySet().join(",")
}

// ArrayList Expando (see processGFF)
ArrayList<Expando> featureList = processGFF(gffFile, ipdProtMap, ipdNucMap,
                                            ipdGeneMap, descSeqMap)
if(debugging <= 3) { 
    err.println "${featureList.size()} features"
}

outFile = new File(gffFile.toString().replaceFirst("_augustus.gff", "_mod.gff"))
gffWriter = new PrintWriter(outFile.newOutputStream(), true)
// output
featureList.each { vList ->
    //(id, gene, seqStartIndex, seqEndIndex) = parseDescription(desc)
    if(debugging <= 1) { 
        err.println "outputting to ${outFile}"
    }
    vList.each { v ->
        if(debugging <= 1) { 
            err.println "output feature ${v}"
        }
        correctAugustus(v, ipdNucMap, ipdGeneMap) // reannotate 2DP1 and 3DP1
        if(debugging <= 2) { 
            err.println "after correctAugustus full match, gene now ${v.gene}, pcall now ${v.pcall}, ccall now ${v.ccall}, gcall now ${v.gcall}"
        }
        outputGFF(v, ipdNucMap, gffWriter)
    }
} // each feature
gffWriter.close()
outputGLs(featureList, outGLFile)

/*
 * outputGFF
 * 
 * @param ipdNucMap Maps of description -> sequence for '_nuc.fasta' file from IPD-KIR
 */
def void outputGFF(Expando query, Map<String,String> ipdNucMap, PrintWriter gbOut) {
    if(debugging <= 1) {
        err.println "outputGFF(${query.gene})"
//        err.println query
    }
    Expando eFormat = initGFF()

    // get best IPD-KIR call
    String bestGene = null
    String bestAllele = null
    String bestRes = null
    gcall = query.gcall ? query.gcall : ""
    ccall = query.ccall ? query.ccall : ""
    pcall = query.pcall ? query.pcall : ""
    if(!gcall.contains("NEW")) {
        bestAllele = gcall
        bestRes = "full gene"
    } else if(!ccall.contains("NEW")) {
        bestAllele = ccall
        bestRes = "cDNA"
    } else if(!pcall.contains("NEW")) {
        bestAllele = pcall
        bestRes = "protein"
    }
    if(bestAllele != null) { 
        (bestGene, allele) = bestAllele.split('\\*')
    }

    ArrayList<String> gffLines = query.newGFF
    partial = true // marking partial if 7 exons or less
    for(int i; i < gffLines.size(); i++) { 
        l = gffLines[i].join('\t')
        if(l.contains("CDS8")) {
            partial = false
            break
        }
    }
    //err.println "partial=${partial}"//remove(todo)
    gffLines.each { l ->
	    if(debugging <= 2) {
		    err.println "outputGFF: l=${l}"
	    }
        desc = l[eFormat.DESC]    // e.g., KP420440_2DL2L3S3S4S5_18152-34641
        type = l[eFormat.TYPE]
        start = l[eFormat.START].toInteger()
        end = l[eFormat.END].toInteger()    
        att = l[eFormat.ATTRIBUTES]

        // type changes
        // replace transcription_start_site and transcription_end_site
        type = type.replaceFirst("transcription_start_site", "five_prime_UTR")
        type = type.replaceFirst("transcription_end_site", "three_prime_UTR")

        // lift the feature coordinates to the full contig/haplotype sequence
        // e.g., KP420440_2DL2L3S3S4S5_18152-34641
        coordStartIndex = desc.lastIndexOf('_')
        (coordStart, coordEnd) = desc[coordStartIndex+1..-1].split('-')
        Integer coordStart = coordStart.toInteger()
        Integer coordEnd = coordEnd.toInteger()
	    if(debugging <= 2) {
		    err.println "outputGFF: desc=${desc}, coordStart=${coordStart}, coordEnd=${coordEnd}"
            err.println "outputGFF: att=${att}"
	    }
        start = start + coordStart
        end = end + coordStart
	    if(debugging <= 2) {
		    err.println "outputGFF: new feature start=${start}, new feature end=${end}"
	    }
        descTmp = desc[0..coordStartIndex-1]
        locusIndex = descTmp.lastIndexOf('_')
        initialDesc = descTmp[0..locusIndex-1]
        descLocus = descTmp[locusIndex+1..-1]

        // add the start location to the gene name to make it unique
        semicolonIndex = att.indexOf(';')
        dotIndex = att.indexOf('.')
        firstIndex = semicolonIndex
        if(dotIndex < semicolonIndex) {
            firstIndex = dotIndex
        }
        id = att[0..firstIndex-1].replaceFirst("ID=", "")
        idNew =  id + "-${coordStart}"
        att = att.replaceAll(id, idNew)
        //err.println "id=${id}, type=${type}, idNew=${idNew}, new att: ${att}"//todo
        // add IPD-KIR gene and allele annotation to gene type
        if(type == "gene") {
            gbOut.println "# Feature:${initialDesc}-${coordStart}"
            if(bestAllele != null) { 
                att += "; gene=${bestGene}; allele=${bestAllele}"
                att += "; note=${KIR_NOMEN_VER}"
            }
            if(query.gene.contains("2DP1") || query.gene.contains("3DP1") ||
               (partial == true)) {
                // http://www.insdc.org/documents/pseudogene-qualifier-vocabulary
                att += "; pseudogene=unprocessed"
            }
            att += "; locus_tag=KIR_${coordStart}"
        } else if(type == "mRNA") {
            att += "; product=Killer cell immunoglobulin-like receptor"
        }

        gbOut.println "" + [initialDesc, l[eFormat.AUGUSTUS], type, start, end,
                            l[eFormat.SCORE], l[eFormat.STRAND], l[eFormat.PHASE],
                            att].join('\t')
    }
    if(debugging <= 1) {
        err.println "outputGFF: return"
    }
    return
} // outputGFF

/*
 * initGFF
 *
 * @return an Expando with the following three properties; the value of each
 * is the Java (0-based) index of its column location in the GFF file: 
 * TYPE, START, END, SCORE, PHASE, ATTRIBUTES.
 * http://gmod.org/wiki/GFF#GFF3_Format
 * 
 */
def Expando initGFF() {
    e = new Expando()
    e.DESC = 0
    e.AUGUSTUS = 1
    e.TYPE = 2
    e.START = 3
    e.END = 4
    e.SCORE = 5
    e.STRAND = 6
    e.PHASE = 7
    e.ATTRIBUTES = 8
    return e
} // initGFF

/*
 * processGFF 
 *
 * Loads the GFF file and puts it into a map (ID -> Expando). The GFF is 0-based
 * indexes. Convert to the 1-based that GenBank wants.
 *
 * @param gffFile the GFF output file as produced by AugustusWeb
 * @param ipdProtMap Maps of description -> sequence for  '_prot.fasta' file from IPD-KIR
 * @param ipdNucMap Maps of description -> sequence for '_nuc.fasta' file from IPD-KIR
 * @param ipdGenMap Maps of description -> sequence for '_gen.fasta' file from IPD-KIR
 * @param descSeqMap Map from description to DNASequence from input fasta file
 * @return a List of Expandos, each containing:
 *   id: id of the individual
 *   gene; the IPD name of the gene
 *   desc: the description of the sequence from the fasta
 *   augustusGeneID: AUGUSTUS' gene and transcript id; e.g., g1.t1
 *   sequence; the full sequence from the fasta as DNASequence
 *   startSeq: 5' start of the input sequence (1-based index)
 *   endSeq: 3' end of the input sequence (1-based index)
 *   start5p: 5' start of the KIR gene (1-based index)
 *   end3p: 3' end of the KIR gene (1-based index)
 *   pcall: the ultimate protein call/interpretation as a gl string (=~ 'NEW' if new)
 *   ccall: the ultimate cDNA call/interpretation as a gl string (=~ 'NEW' if new)
 *   gcall: the ultimate dna full gene call/interpretation as a gl string (=~ 'NEW' if new)
 *   exonStartTreeSet = [:] // exon index -> start of the exon (1-based index)
 *   exonEndTreeSet = [:] // feature -> end of the exon (1-based index)
 *   newGFF = ArrayList<String> of the new GFF lines
 */
ArrayList<Expando> processGFF(File gffFile, Map<String,String> ipdProtMap,
                              Map<String,String> ipdNucMap,
                              Map<String,String> ipdGeneMap,
                              LinkedHashMap<String, DNASequence> descSeqMap) {
    if(debugging > 3) {
        err.println "processGFF(gffFile=${gffFile})"
    }

    Expando eFormat = initGFF() // properties are column properties

    // this will be returned; list of Expando
    ArrayList<Expando> geneList = new ArrayList()
    ICsvListReader listReader = null
    listReader = new CsvListReader(new FileReader(gffFile),
                                   CsvPreference.TAB_PREFERENCE)
    // load the Augustus output fasta files for the amino acid (.aa),
    // exon dna (.cdsexons), and coding seq (.codingseq)
    LinkedHashMap<String, ProteinSequence> gffAAMap
    LinkedHashMap<String, DNASequence> gffCDSMap
    LinkedHashMap<String, DNASequence> gffCodingMap
    (gffAAMap, gffCDSMap, gffCodingMap) = loadGFFMaps(gffFile)

    List<Object> line
    while( (line = listReader.read()) != null ) {  // reading GFF
        if(line.size() == 1) { // ignore comments
            continue;
        }
        /* could put phase processing here
         err.println "phase ${line[eFormat.PHASE]}"
         if(line[eFormat.PHASE] == "-") { // skip reverse phase for now
         err.println "skipping phase ${line[eFormat.PHASE]}"
         continue
     }
         */
        // TYPE, START, END, SCORE, PHASE, ATTRIBUTES
        String type = line[eFormat.TYPE]
        if(debugging <= 1) { 
            err.println "processGFF: line=${line}"
        }
        // if the gene isn't all there, transcription_start_site will be absent
        // and the order will be gene, mRNA, intron, CDS, exon...
        if(type != "gene") {
            continue;
        }
        while((line != null) && (line != "")) {
            Expando eGene = null
            (eGene, line) = processGFFGene(eFormat, line, listReader, descSeqMap, gffAAMap,
                                           gffCDSMap, gffCodingMap, ipdProtMap, ipdNucMap, ipdGeneMap)
            if((eGene == null) || eGene.newGFF.size() == 0) {
                break
            }
            previousEntryList = geneList
            ArrayList<Expando> toAdd = new ArrayList()
            ArrayList<Expando> toRemove = new ArrayList()
            if(previousEntryList.size() > 0) { 
                previousEntryList.each { previousEntry ->
                    // the second+ mRNA won't have a gene line
                    if(eGene.newGFF[0][eFormat.TYPE] != "gene") {
                        // add the parent from the last entry
                        lastEntry = previousEntryList[-1]
                        // add the gene line to the new
                        eGene.newGFF.add(0, lastEntry.newGFF[0])
                    }
                    if(debugging <= 1) {
                        if(previousEntry != null) { 
                            err.println "processGFF: previousEntry=${previousEntry}"
                            err.println "processGFF: previous ccall=${previousEntry.ccall ? previousEntry.ccall : "null"}"
                        }
                        err.println "eGene ccall2=${eGene.ccall}"
                    }
                    // replace NEWs with something that follows
                    // todo: improve this?; use scores?
                    /* only needed if multiple transcripts are evaluated?
                     if((previousEntry == null) || (previousEntry.gcall == NEW_ALLELE)) { 
                     if(debugging <= 2) {
                     err.println "processGFF: adding ${eGene.gene}"
                 }
                     geneMap[eGene.id] = eGene//todo
                 }*/
                    if(GENE_SYSTEM == "KIR") {
                        (toAdd, toRemove) = newKIREntryIsBetter(previousEntry, eGene, e.exonStartTreeSet, e.exonEndTreeSet)
                    } // if KIR
                    if(debugging <= 1) { 
                        err.println "processGFF ${eGene.gene}=${eGene.ccall}"
                    }
                }
                if(debugging <= 2) { 
                    err.println "${eGene.id} removing " + toRemove
                    err.println "${eGene.id} adding " + toAdd
                }
                previousEntryList.removeAll(toRemove)
                previousEntryList.addAll(toAdd)
            } else {
                if(debugging <= 2) { 
                    err.println "processGFF: adding ${eGene.id} " + eGene
                }
                geneList.add(eGene)
            }
        } // repeat processGFFGene until it returns blank or null line
    } // each row

    if(debugging <= 3) { 
        err.println "processGFF: return ${geneList.size()} IDs"
    }
    return geneList
} // processGFF

/*
 * Parses a single gene section of a GFF file.
 *
 * @param eFormat an Expando; see initGFF()
 * @param line ArrayList from a CsvListReader
 * @param listReader a CsvListReader containing the remaining lines from the GFF file
 * @param ipdDir File of the directory of the '_prot.fasta', '_nuc.fasta',
 *             and '_gen.fasta' files from IPD-KIR
 * @param descSeqMap Map from description to DNASequence from input fasta file
 * @param gffAAMap Map of GFF ProteinSequences
 * @param gffCDSMap Map of GFF exonic DNASequences
 * @param gffCodingMap Map of GFF full cDNA DNASequences
 * @return an Expando (see processGFF) and the current line
 * @todo modularize this (see exon scope)
 */
def ArrayList<Object> processGFFGene(Expando eFormat, ArrayList line,
                                     CsvListReader listReader, LinkedHashMap<String, DNASequence> descSeqMap,
                                     LinkedHashMap<String, ProteinSequence> gffAAMap,
                                     LinkedHashMap<String, DNASequence> gffCDSMap,
                                     LinkedHashMap<String, DNASequence> gffCodingMap,
                                     Map<String,String> ipdProtMap,
                                     Map<String,String> ipdNucMap,
                                     Map<String,String> ipdGeneMap) {
    if(debugging <= 1) { 
        err.println "processGFFGene(eFormat=${eFormat}, line=${line})"
    }
    e = new Expando() // return value
    e.exonStartTreeSet = new TreeSet()
    e.exonEndTreeSet = new TreeSet()
    e.newGFF = new ArrayList()
    if((line[eFormat.TYPE] != null) && (line[eFormat.TYPE] == "gene")) { 
        e.newGFF.add(line)
        line = listReader.read() // get past the "gene" line
    }
    if(debugging <= 1) { 
        err.println "processGFFGene: new line=${line}"
    }
    HashSet<String> pCallSet = new HashSet()
    HashSet<String> cCallSet = new HashSet()
    HashSet<String> gCallSet = new HashSet()
    String gene // GFF gene name (e.g., 'g1')
    String transcript // GFF transcript name (e.g., 't1')
    Integer exonIndex = 0 // 0-based Array indexes
    while(line != null) {
        if(line.size() == 1) { // ignore comments
            line = listReader.read()
            continue;
        }

        String type = line[eFormat.TYPE]
        if(debugging <= 1) {
            err.println "processGFFGene: type=${type}"
        }
        Integer start = line[eFormat.START].toInteger()
        Integer end = line[eFormat.END].toInteger()
        if((type == "transcription_start_site") ||
           ((transcript == null) && (type == "mRNA"))) {
            // get the ID, and location of the sequence within the haplotype
            // e.g., cB02~tB01_GU182355.1_2DL4_67155-80842
            e.start5p = start
            desc = line[eFormat.DESC]
            e.desc = desc
            (id, gene, seqStartIndex, seqEndIndex) = parseDescription(desc)
            e.startSeq = seqStartIndex
            e.endSeq = seqEndIndex
            e.id = id
            faSeq = descSeqMap[e.desc]
            if(debugging <= 2) {
                err.println "processGFFGene: setting sequence for ID ${e.id}, desc=${desc}"
            }
            if(faSeq == null) {
                err.println "ERROR: couldn't find sequence in fasta for description ${e.desc}"
                err.println descSeqMap.keySet().sort().join(",")
                System.exit(1)
            }
            e.sequence = faSeq // DNASequence
            
            attributes = line[eFormat.ATTRIBUTES]
            if(debugging <= 2) {
                err.println "processGFFGene: attributes=${attributes}"
            }
            // e.g., ID=g1.t1.tss1;Parent=g1.t1;
            (prefixAtt, rest) = attributes.split(';')
            (gene, transcript) = prefixAtt.split('\\.')
            gene = gene.replaceFirst("ID=", "")
            ipdGene = GENE_NAME
            if(debugging <= 3) {
                err.println "processGFFGene: ${gene}.${transcript}"
                err.println "processGFFGene: ipdGene=${ipdGene}"
            }
            e.augustusGeneID="${gene}.${transcript}"   // e.g., g1.t1
            e.gene = GENE_NAME
            e.newGFF.add(line)
        } else if(type == "CDS") {
            e.exonStartTreeSet.add(start)
            e.exonEndTreeSet.add(end)
            if(debugging <= 3) {
                err.println "processGFFGene: gene ${gene} exon ${exonIndex}"
            }
            // this is what happens when the 5' part of the gene is missing
            // todo: handle this better
            if(gene == null) {
                if(debugging <= 1) {
                    err.println "processGFFGene: breaking"
                }
                break;
            }
            if(exonIndex == 0) {
                // process all the exons
                e = processGFFExon(gene, transcript, ipdProtMap, ipdNucMap,
                                   ipdGeneMap, gffAAMap, gffCDSMap, gffCodingMap, e,
                                   pCallSet, cCallSet)
                e = processGFFFullGene(gene, ipdGeneMap, e, gCallSet, cCallSet, pCallSet)
                if(debugging <= 1) { 
                    err.println "processGFFGene: pcall set, size now2 ${pCallSet.size()}"
                    err.println "processGFFGene: ccall set, size now2 ${cCallSet.size()}"
                }
            }
            exonIndex++;
            e.newGFF.add(line)
        } else if(type == "mRNA") {
                e.newGFF.add(line)
        } else if((type == "transcription_end_site") || (type == "gene")) {
            e.end3p = end
            if(debugging <= 2) {
                err.println "processGFFGene: found tts"
            }
            if(type != "gene") {
                e.newGFF.add(line)
                line = listReader.read()
            }                    
            break; // time to returnfrom the method
        }
        // noop; see break at "tts"
        if(type != "gene") {
            line = listReader.read()
        }
    } // each line
    if(debugging <= 3) {
        err.println "processGFFGene: end of gff pCallSet size=${pCallSet.size()}, cCallSet size=${cCallSet.size()}"
    }            
    if(pCallSet.isEmpty()) {
        e.pcall = "${GENE_NAME}*${NEW_ALLELE}"
    } else { 
        e.pcall = pCallSet.sort().join("/")
    }
    if(cCallSet.isEmpty()) {
        e.ccall = "${GENE_NAME}*${NEW_ALLELE}"
    } else { 
        e.ccall = cCallSet.sort().join("/")
    }
    if(gCallSet.isEmpty()) {
        e.gcall = "${GENE_NAME}*${NEW_ALLELE}"
    } else { 
        e.gcall = gCallSet.sort().join("/")
    }
    
    if(debugging <= 3) {
        err.println "processGFFGene: return ${e.pcall}, ${e.ccall}, ${e.gcall}, line=${line}"
    }
        return [e, line]
} // processGFFGene

/*
 * Parses and interprets a exon section of a GFF file.
 *
 * @param gene the GFF gene name (e.g., 'g1')
 * @param transcript the GFF transcript number (e.g., 't1')
 * @param ipdProtMap the IPD _prot.fasta sequences
 * @param ipdNucMap the IPD _nuc.fasta sequences
 * @param ipdGeneMap the IPD _gen.fasta sequences
 * @param gffAAMap Map of GFF ProteinSequences
 * @param gffCDSMap Map of GFF exonic DNASequences
 * @param gffCodingMap Map of GFF full cDNA DNASequences
 * @param e to fill in; see processGFF
 * @param pCallSet the protein call(s) for this gene
 * @param cCallSet the cDNA call(s) for this gene
 * @return the Expando that was input; see processGFF
 */
def Expando processGFFExon(String gene, String transcript,
                           Map<String, String> ipdProtMap, 
                           Map<String, String> ipdNucMap, 
                           Map<String, String> ipdGeneMap,
                           LinkedHashMap<String, ProteinSequence> gffAAMap,
                           LinkedHashMap<String, DNASequence> gffCDSMap,
                           LinkedHashMap<String, DNASequence> gffCodingMap,
                           Expando e, HashSet<String> pCallSet,
                           HashSet<String> cCallSet) {
    if(debugging <= 1) {
        err.println "processGFFExon(gene=${gene}, transcript=${transcript})"
    }
    // check protein sequence (AUGUSTUS: .aa & IPD: _prot.fasta)
    m = "${gene}.${transcript}"
    gffSeq = gffAAMap[m]
    if((gffSeq == null) || (gffSeq == "")) {
        err.println "ERROR: processGFFExon: couldn't find ${m} in gffAAMap for ${e.id}"
        err.println "processGFFExon: gffAAMap keys " + gffAAMap.keySet().join(",")
    }
    if(debugging <= 1) {
        err.println "processGFFExon: seq ${m}=${gffSeq}"
    }

    ipdProtMap.each { header, seq ->
        if(debugging <= 1) {
            err.println "processGFFExon: protein ${header}=${seq}"
        }
        protTrue = sequenceEquals(gffSeq, seq)
        ipdGene = GENE_NAME
        // protein match *
        if(protTrue) {
            // check cDNA sequence (AUGUSTUS: .codingseq & IPD: _nuc.fasta)
            ipdGeneNew = ipdGene
            if(debugging <= 2) { 
                err.println "processGFFExon: ipdGene 1 =${ipdGene}"
            }
            (ipdGeneNew, pcall) = truncateToProtein(header)
            e.gene = ipdGeneNew
            fullName = "${ipdGeneNew}*${pcall}"
            if(debugging <= 3) {
                err.println "processGFFExon: *protein match: ${fullName}"
            }
            if(debugging <= 1) {
                err.println "processGFFExon: gffSeq=${gffSeq}"
                err.println "processGFFExon: ipd seq=${seq}"
            }
            pCallSet.add(fullName)
            // find the AUGUSTUS cDNA sequence
            gffCodingMap.each { cHeader, cSeq ->
/*                if(debugging <= 2) {
                    err.println "processGFFExon: ${cHeader}"
                }*/
                if(cHeader.endsWith(m)) {
                    (ipdGeneRet, ccall, fullName) = getcDNAMatches(cSeq, ipdNucMap)
                    if(ipdGeneRet != null) {
                        e.gene = ipdGeneRet
                        cCallSet.add(fullName)
                        if(debugging <= 1) { 
                            err.println "processGFFExon: adding ${fullName} to ${e.id} ccall set, size now ${cCallSet.size()}"
                        }
                    }
                } // found the correct sequence in the .codingseq file
            } // each sequence from the AUGUSTUS .codingseq file
        } // protein match
    } // each ipd protein sequence for this gene

    if(debugging <= 1) {
        err.println "processGFFExon: return ${e.pcall}, ${e.ccall}, pCallSet size=${pCallSet.size()}, cCallSet size=${cCallSet.size()}"
    }
    return e
} // processGFFExon

/*
 * Annotates the 'full gene' sequence wrt IPR alleles.
 *
 * @param gene the GFF gene name (e.g., 'g1')
 * @param ipdGeneMap the IPD _gen.fasta sequences
 * @param e to fill in; see processGFF
 * @param gCallSet the full-sequence allele call(s) for this gene
 * @param cCallSet the cDNA call(s) for this gene
 * @param pCallSet the protein call(s) for this gene
 * @return the Expando that was input; see processGFF
 */
def Expando processGFFFullGene(String gene, Map<String, String> ipdGeneMap,
                               Expando e, Set gCallSet, Set cCallSet, Set pCallSet) {
    if(debugging <= 1) {
        err.println "processGFFFullGene(id=${e.id}, gene=${gene})"
    }
    if(e.sequence == null) {
        err.println "processGFFFullGene: ERROR: no sequence for id=${e.id}, gene=${gene}"
        System.exit(1)
    }
    genSeq = e.sequence.getSubSequence(e.start5p, e.end3p).getViewedSequence()
    genSeqStr = genSeq.getSequenceAsString()
    if(debugging <= 3) {        
        err.println "processGFFFullGene: ${gene} sequence is " +
            "${genSeq.getLength()} bp long"
    }
    
    ipdGene = GENE_NAME
    ipdGeneMap.each { header, seq ->
        // full sequence match *
        if(sequenceContains(seq, genSeq)) {
            (ipdGeneNew, gcall) = interpFullGene(header)
            fullName = "${ipdGeneNew}*${gcall}"
            if(debugging <= 3) {
                err.println "processGFFFullGene: *match on full gene: ${fullName}"
            }
            gCallSet.add(fullName)
            (cGene, ccall) = truncateToCDNA(fullName)
            fullCDNA = "${cGene}*${ccall}"
            cCallSet.add(fullCDNA)
            (ipdGeneNew, pcall) = truncateToProtein(fullCDNA)
            fullProt = "${ipdGeneNew}*${pcall}"
            pCallSet.add(fullProt)
            e.gene = ipdGeneNew
        }
    } // each ipd/imgt full sequence

    if(!gCallSet.isEmpty()) {
        e.gcall = gCallSet.sort().join("/")
    } else {
        e.gcall = "${ipdGene}*${NEW_ALLELE}"
    }
    
    if(debugging <= 1) {
        err.println "processGFFFullGene: return ${e.gcall}"
    }
    return e
} // processGFFFullGene

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


/*
 * Loads the three relevant fasta files from the AUGUSTUSWweb output.
 * Amino acid (.aa), exon dna (.cdsexons), and coding seq (.codingseq).
 *
 * @param gffFile the GFF output file as produced by AugustusWeb
 * @return a List containing three LinkedHashMaps:
 *                    0: <String, ProteinSequence> for .aa
 *                    1: <String, DNASequence> for .cdsexons
 *                    2: <String, DNASequence> for .codingseq
 */
def List<LinkedHashMap> loadGFFMaps(gffFile) { 
    LinkedHashMap<String, ProteinSequence> gffAAMap
    LinkedHashMap<String, DNASequence> gffCDSMap
    LinkedHashMap<String, DNASequence> gffCodingMap

    // .aa
    protFileName = gffFile.getAbsolutePath().replace(".gff", ".aa")
    if(debugging <= 3) { 
        err.println "loadGFFMaps: loading ${protFileName}..."
    }
    gffAAMap = new LinkedHashMap()
    pf = new File(protFileName)
    if(pf.length() != 0) { 
        gffAAMap = FastaReaderHelper.readFastaProteinSequence(pf)
    }
    if(debugging <= 3) { 
        err.println gffAAMap.keySet().size() + " items"
    }

    // .cdsexons
    nucFileName = gffFile.getAbsolutePath().replace(".gff", ".cdsexons")
    if(debugging <= 3) { 
        err.println "loadGFFMaps: loading ${nucFileName}..."
    }
    gffCDSMap = new LinkedHashMap()
    nf = new File(nucFileName)
    if(nf.length() != 0) { 
        gffCDSMap = FastaReaderHelper.readFastaDNASequence(nf)
    }
    if(debugging <= 3) { 
        err.println gffCDSMap.keySet().size() + " items"
    }
    
    // .codingseq
    nucFileName = gffFile.getAbsolutePath().replace(".gff", ".codingseq")
    if(debugging <= 3) { 
        err.println "loadGFFMaps: loading ${nucFileName}..."
    }
    gffCodingMap = new LinkedHashMap()
    nf = new File(nucFileName)
    if(nf.length() != 0) { 
        gffCodingMap = FastaReaderHelper.readFastaDNASequence(nf)
    }
    if(debugging <= 3) { 
        err.println gffCodingMap.keySet().size() + " items"
    }
    
    return [gffAAMap, gffCDSMap, gffCodingMap]
} // loadGFFMaps

/* 
 * parseDescription 
 *
 * Gets the ID, and location of the sequence within the haplotype. The format
 * is nomenclature_ID_gene_start-end. Assuming 0-based index
 * e.g., cB02~tB01_GU182355.1_2DL4_67155-80842
 * 
 * @param String the description from the fasta files
 * @return ArrayList of nomenclature, ID (String), sequence start index (Integer), end index (Integer); or null values in the Array
 */
def ArrayList parseDescription(String desc) {
    if(debugging <= 2) {
        err.println "parseDescription(desc=${desc})"
    }
    ArrayList retArray = new ArrayList(3)

    // nomenclature_ID_gene_start-end
    // e.g., cB02~tB01_GU182355.1_2DL4_67155-80842
    // struct _ ID _ gene _ start-end
    Queue<String> descQ = new LinkedList() as Queue
    desc.split('_').each { descQ.push(it) }
    locations = descQ.pop()
    gene = descQ.pop()
    id = descQ.join('_') + "_${gene}_${locations}"
    if(debugging <= 2) {
        //err.println "matcher " + matcher 
        err.println "id " + id 
        err.println "gene " + gene 
        err.println "locations " + locations 
    }

    (sString, eString) = locations.split("-")
    Integer start = sString.toInteger() 
    Integer end = eString.toInteger()
    
    if(debugging <= 2) {
        err.println "parseDescription: id=${id}, gene=${gene}, " +
            "gene=${gene}, start=${start}, end=${end}"
    }

    return [id, gene, start, end]
} // parseDescription

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
def ArrayList<String> interpFullGene(String allele) {
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
 * Outputs the GL strings to a file. One for cDNA resolution & one for full.
 *
 * @param featureList a List of Expandos (see processGFF())
 * @param outGLFile the file to which to write the GL strings
 * @todo rename: not gl any more
 */
def void outputGLs(ArrayList<Expando> featureList, File outGLFile) {
    if(debugging <= 1) {
        err.println "outputGLs(outGLFile=${outGLFile.getName()})"
    }
    FileWriter outGLWriter = new FileWriter(outGLFile)

    outGLWriter.println "ID\taugID\tgene\tprotein\tcDNA\tfull"
    featureList.each { e ->
        outGLWriter.println "${e.id}\t${e.augustusGeneID}\t${e.gene}\t${e.pcall}\t${e.ccall}\t${e.gcall}"
    }
        
    outGLWriter.close()
    if(debugging <= 1) {
        err.println "outputGLs: return"
    }    
} // outputGLs


/*
 * handleArgs
 * 
 * @param args the command line arguments
 * @return list of command line arguments in a usable form; use -h for more info
 */
def ArrayList handleArgs(String[] args) {
    USAGE = "gff2ftGene.groovy [args]"
    Switch help = new Switch("h", "help", "display help message")
    ipdArg = new FileArgument("i", "ipd",
                              "location of IPD-KIR directory", true)
    gffArg = new FileArgument("f", "gff", "GFF AUGUSTUS output file",
                              true)
    seqArg = new FileArgument("s", "seq", "fasta file containing the sequence",
                              true)
    outGLArg = new FileArgument("l", "gl", "output gl file", true)
    geneSysArg = new StringArgument("g", "gene", "e.g., 'KIR2DL4' or 'HLA-A'",
                                    true)
    ArgumentList arguments = new ArgumentList(help, ipdArg, gffArg, seqArg,
                                              outGLArg, geneSysArg)

    CommandLine commandLine = new CommandLine(args)
    try {
        CommandLineParser.parse(commandLine, arguments)
        if (help.wasFound()) {
            Usage.usage(USAGE, null, commandLine, arguments, System.out)
            System.exit(-2)
        }
    } catch (CommandLineParseException e) {
        Usage.usage(USAGE, null, commandLine, arguments, System.err)
        System.exit(-1)
    }

    return [ipdArg.getValue(), gffArg.getValue(), seqArg.getValue(),
            outGLArg.getValue(), geneSysArg.getValue()]
} // handleArgs

// todo: document
def Boolean sequenceEquals(org.biojava.nbio.core.sequence.template.Sequence a,
                           org.biojava.nbio.core.sequence.template.Sequence b) {
    if(debugging <= 1) {
        err.println "sequenceEquals()"
    }
    Boolean ret = false
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
 * Returns new is the transcript in eGene is a better call than
 * previousEntry. This is currently done by roughly checking the
 * expected sizes of the exons.
 */
def List<List<Expando>> newKIREntryIsBetter(Expando previousEntry, Expando eGene, TreeSet<String> exonStartTreeSet,
                                TreeSet<String> exonEndTreeSet) {
    if(debugging <= 1) { 
        err.println "newKIREntryIsBetter()"
    }
    Boolean ret = false

    List<Expando> toAdd = new ArrayList()
    List<Expando> toRemove = new ArrayList()
    int previousSize = previousEntry.exonEndTreeSet.size()
    int newSize = eGene.exonEndTreeSet.size()
    pid = previousEntry.id
    gid = eGene.id
    pgene = previousEntry.gene
    ggene = eGene.gene
    Integer prevStart = previousEntry.startSeq
    Integer geneStart = eGene.startSeq
    Integer prevEnd = previousEntry.endSeq
    Integer geneEnd = eGene.endSeq
    if(debugging <= 1) { 
        err.println "newKIREntryIsBetter: previous id=${pid}, eGene id=${gid}, previous gene=${pgene}, eGene gene=${ggene}, previousSize=${previousSize}, newSize=${newSize}"
        err.println "newKIREntryIsBetter: prevStart=${prevStart}, prevEnd=${prevEnd}, geneStart=${geneStart}, geneEnd=${geneEnd}"
    }

    // check if it is a different gene first
    //old if((pid != gid) || (pgene != ggene) || ((prevStart >= geneStart) && (prevStart <= geneEnd)) || ((geneStart >= prevStart) && (geneStart <= prevEnd)) ) {
    if((pid == gid) &&
       ((prevStart >= geneStart) && (prevStart <= geneEnd)) || ((geneStart >= prevStart) && (geneStart <= prevEnd)) ) {
        if(debugging <= 3) {
            err.println ((prevStart >= geneStart) && (prevStart <= geneEnd))
            err.println ((geneStart >= prevStart) && (geneStart <= prevEnd))
            err.println "newKIREntryIsBetter: same id/gene pair: previous id=${pid}, eGene id=${gid}, previous gene=${pgene}, eGene gene=${ggene}"
        }
        // prioritize the gene calls
        Iterator siter = eGene.exonStartTreeSet.iterator()
        Iterator eiter = eGene.exonEndTreeSet.iterator()
        Iterator psiter = previousEntry.exonStartTreeSet.iterator()
        Iterator peiter = previousEntry.exonEndTreeSet.iterator()
        for(; siter.hasNext() && psiter.hasNext();) {
            Integer start = siter.next()
            Integer end = eiter.next()
            Integer pstart = psiter.next()
            Integer pend = peiter.next()
            Integer previousExonLength = pend - pstart + 1
            Integer exonLength = end - start + 1

            previousContains = KIR_EXON_SIZES.contains(previousExonLength)
            contains = KIR_EXON_SIZES.contains(exonLength)

            if(!previousContains && contains) {
                toAdd.add(eGene)
                toRemove.add(previousEntry)
            }
        } // each exon
        if(debugging <= 1) { 
            err.println "newKIREntryIsBetter: return; toAdd=${toAdd}, toRemove=${toRemove}"
        }
        return [toAdd, toRemove]
    } else {
        toAdd.add(eGene)
    }
    
    if(debugging <= 1) { 
        err.println "newKIREntryIsBetter: return; toAdd=${toAdd}, toRemove=${toRemove}"
    }
    return [toAdd, toRemove]
} // newKIREntryIsBetter

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
 * Adjust the AUGUSTUS annotation for the two KIR pseudo genes, 2DP1 and 3DP1.
 * Change ccall and pcall in the expando.
 * 
 * @param v an Expando (see ProcessGFF)
 * @param ipdNucMap Maps of description -> sequence for '_nuc.fasta' file from IPD-KIR
 * @todo: expand to full sequence
 * @todo: rename: not just for pseudo
 */ 
def void correctAugustus(Expando v, Map<String, String> ipdNucMap, Map<String, String> ipdGeneMap) {
    gene = v.gene
    if(debugging <= 1) {
        err.println "correctAugustus(ID=${v.id}, gene=${v.gene})"
    }
    if(!gene.contains("2DP1") && !gene.contains("3DP1")) {
        return
    }

    cDNANew = new String() // the new cDNA sequence
    numExons = v.exonStartTreeSet.size()//left off: wrong here?
    if(debugging <= 2) {
        err.println "correctAugustus: ${numExons} exons"
    }
    // record the changes (indexes) to be made to exonStartTreeSet and exonEndTreeSet
    ArrayList<Integer> startToAdd = new ArrayList()
    ArrayList<Integer> startToRemove = new ArrayList()
    ArrayList<Integer> endToAdd = new ArrayList()
    ArrayList<Integer> endToRemove = new ArrayList()
    Iterator siter = v.exonStartTreeSet.iterator()
    Iterator eiter = v.exonEndTreeSet.iterator()    
    for(i=0; siter.hasNext(); i++) { 
        Integer start = siter.next()
        Integer end = eiter.next()
        Integer length = end - start + 1
        newStart = start
        newEnd = end
        if(debugging <= 2) {
            err.println "correctAugustus: exon ${i}, start=${start}, end=${end}, length=${length}"
        }

        /*remove(todo)
        // skip an extra exon that augustus calls  2DP1 (241 bp)
        if(v.gene.contains("2DP1") && (length == 241)) {//todo: move this
            startToRemove.add(start)
            endToRemove.add(end)
            if(debugging <= 2) {
                err.println "correctAugustus: exon ${i}, removing start=${start}, end=${end}, length=${length}"
            }
        }*/
        if((i+1) == numExons) { // for last exon, include 3' UTR
            // augustus has the correct end, but not start of last exon for 2DP1
            // augustus has the correct start, but not end of last exon for 3DP1
            if(v.gene.contains("2DP1")) {
                newStart = end - 177 + 1
                // could just add ~264 to start index based on manual observation
                seq = "GTTGTCTCCTGCCCATGA"
                end3p = v.sequence.getSequenceAsString().indexOf(seq)
                if(end3p < 0) {
                    err.println "correctAugustus: WARNING: couldn't find 3' UTR in 2DP1 for ${v.id}"
                } else { 
                    newEnd = end3p + seq.length()
                    if(debugging <= 2) {
                        err.println "correctAugustus: 2DP1 3' UTR=${end3p}"
                    }
                    v.end3p = newEnd
                }

                //newEnd = v.end3p
                //todo(keep?)newEnd = newStart + 177 - 1
            } else if(v.gene.contains("3DP1")) { 
                newEnd = start + 294 - 1
            }
            if(debugging <= 2) {
                err.println "correctAugustus: removing ${start}, ${end}"
                err.println "correctAugustus: adding ${newStart}, ${newEnd}"
            }
            if(start != newStart) { 
                startToRemove.add(start)
                startToAdd.add(newStart)
            }
            if(end != newEnd) { 
                endToRemove.add(end)
                endToAdd.add(newEnd)
            }
        } else if(i == 0) { // for exon 1, include 5' UTR
//            newStart = v.start5p
            if((newEnd - newStart + 1) == 13) {
                newStart = newEnd - 33
            }
            startToRemove.add(start)
            startToAdd.add(newStart)
        } else if(!KIR_EXON_SIZES.contains(length)) {
            if(debugging <= 4) { 
                err.println "correctAugustus WARNING: skipping exon ${i} length ${length} (${end}-${start}) for ${v.gene}, ${v.id}"
            }
            startToRemove.add(start)
            endToRemove.add(end)
            continue
        }
    } // each exon
    // also, add exon 8 for 2DP1
    if(gene.contains("2DP1")) {
        // exon '8' (including the pseudo exon)
        // the exon is ATGCTGCTGTAATGGACCAAGAGCCTGCAGGGAACAGAACAGCGAATAGCTAG
        // but the 3rd from the last character is variable
        exon8 = "ATGCTGCTGTAATGGACCAAGAGCCTGCAGGGAACAGAACAGCGAATAGC"
        //err.println "correctAugustus: " + v.sequence //todo
        start8 = v.sequence.getSequenceAsString().indexOf(exon8)
        if(start8 < 0) {
            err.println "correctAugustus: WARNING: couldn't find exon 8 in 2DP1 for ${v.id}"
        } else { 
            start8++
            end8 = start8 + 52
            length8 = end8 - start8 + 1
            if(debugging <= 2) {
                err.println "correctAugustus: adding exon 8 ${start8}\t${end8}\t${length8} bp"
            }
            startToAdd.add(start8)
            endToAdd.add(end8)
        }
    } // done adding exon 8

    // remove must come before adding; some may be in both lists
    v.exonStartTreeSet.removeAll(startToRemove)
    v.exonStartTreeSet.addAll(startToAdd)
    v.exonEndTreeSet.removeAll(endToRemove)
    v.exonEndTreeSet.addAll(endToAdd)

    // generate the cDNA sequence from the new exons
    Iterator siter3 = v.exonStartTreeSet.iterator()
    Iterator eiter3 = v.exonEndTreeSet.iterator()
    //err.println v.exonStartTreeSet //todo
    //err.println v.exonEndTreeSet //todo
    for(i=0; siter3.hasNext(); i++) { 
        Integer start = siter3.next()
        Integer end = eiter3.next()
        l3 = end - start + 1
        if(debugging <= 2) {
            err.println "correctAugustus: exon ${i} ${start}\t${end}\t${l3} bp"
        }
        cDNANew += v.sequence.getSequenceAsString()[start-1..end-1]
    } // each exon
    
    // annotate the sequence wrt IPD-KIR
    // annotate wrt cDNA
    cDNANewSeq = new DNASequence(cDNANew, AmbiguityDNACompoundSet.getDNACompoundSet())
    (ipdGeneName, ccall, fullName) = getcDNAMatches(cDNANewSeq, ipdNucMap)
    if(ipdGeneName != null) {
        v.gene = ipdGeneName
        v.ccall = fullName
        (ipdGeneNew, pcall) = truncateToProtein(fullName)
        v.pcall = "${ipdGeneNew}*${pcall}"
        if(debugging <= 2) {
            err.println "correctAugustus: after cDNA match, gene now ${v.gene}, pcall now ${v.pcall}, ccall now ${v.ccall}"
        }
    }
    if(debugging <= 2) {
        err.println "correctAugustus: before checking full gene, v=${v}"
    }
    if(v.end3p == null) { // todo: why is this null here? KP420446 3DP1
        v.end3p = v.sequence.getLength()-1
    }
    // annotate wrt full gene
    ipdGeneMap.each { header, seq ->
        // full sequence match *
        geneSeq = new DNASequence(v.sequence.getSequenceAsString()[v.start5p-1 .. v.end3p-1], AmbiguityDNACompoundSet.getDNACompoundSet())
        if(sequenceContains(seq, geneSeq)) {
            (ipdGeneNew, gcall) = interpFullGene(header)
            v.gene = ipdGeneNew
            fullName = "${ipdGeneNew}*${gcall}"
            if(debugging <= 3) {
                err.println "correctAugustus: *match on full gene: ${fullName}"
            }
            v.gcall = fullName
            (ipdGeneNew, ccall) = truncateToCDNA(header)
            fullName = "${ipdGeneNew}*${ccall}"
            v.ccall = fullName
            (ipdGeneNew, pcall) = truncateToProtein(fullName)
            fullName = "${ipdGeneNew}*${pcall}"
            v.pcall = fullName
            if(debugging <= 2) {
                err.println "correctAugustus: after full match, gene now ${v.gene}, pcall now ${v.pcall}, ccall now ${v.ccall}, gcall now ${v.gcall}"
            }
        }
    } // each ipd/imgt full sequence

} // correctAugustus
