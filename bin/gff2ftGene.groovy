#!/usr/bin/env groovy

/*
 * Convert multi-gene GFF format (as output from chrishah/premaker-plus AUGUSTUS) to 
 *   1. A Feature Table as required for sequence submission to GenBank.
 *      1-based indexes
 *   2. Three GL strings containg IMGT/IPD interpretations: 
 *       protein, cDNA, full gene.
 *
 * e.g., gff2ftGene.groovy -g KIR2DL4 -i ~/doc/KIRDB/2.8.0/fasta/ -f all_2DL4/augustus.gff -s all_seqs_v3_2DL4_features_short.fasta -o seqs_v3_all_2DL4.ft.txt -l seqs_v3_all_2DL4.gl.txt 2> seqs_v3_all_2DL4.ft_err.txt

 * main arguments:
 *   -h  display help message [optional]
 *   -i  location of directory containing IPD-KIR fasta files
 *   -g  gene name (e.g., KIR2DL4 or HLA-A) [required]
 *   -f  location to GFF file output by AUGUSTUS (Web); 0-based indexes
 *   -s, --seq [class java.io.File]  fasta file containing the sequence [required]
 *   -o  location to output GenBank feature table
 *   -l  location to output the GL string
 *
 * Requires 
 *   - Commandline: http://www.dishevelled.org/commandline/
 *   - Super CSV: http://super-csv.github.io/super-csv/
 *   - BioJava 4's Core jar: http://www.biojava/org
 *   - slf4j-api & slf4j-log4j12 http://www.slf4j.org/
 *   - log4j 1 https://logging.apache.org/log4j/1.2/
 *   e.g., 
 *     export CLASSPATH=$HOME/bin/jars/dsh-commandline-1.1.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/super-csv.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/biojava4-core.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/slf4j-api-1.7.5.jar:$CLASSPATH
 *     export CLASSPATH=$HOME/bin/jars/log4j-1.2.15.jar:$CLASSPATH
 *
 * @author Dave Roe
 * @todo indexOf()s need reverse complement checked
 * @todo add pseudogene notes etc
 * @todo replace Commandline ?; see annotateMarkup
 * @todo cB02~tA01_AFA12 BCDcLMHCIJKLFGHCIJKLSCTUVWXYKLSCANJKLFZZZZZGHCIJKLSCANJ 3DL3~2DL1/2DS1/2DS2/2DS2\2DS3~2DL2/2DL3~3DP1~3DP1-2DL4~2DL4~3DL1L2S1~2DS3/2DS5/2DS4~3DL1L2S1
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
(ipdDir, gffFile, seqFile, outFile, outGLFile, GENE_NAME) = handleArgs(args)

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

// ID - ArrayList Expando (see processGFF)
Map<String, ArrayList<Expando>> featureMap = processGFF(gffFile, ipdProtMap, ipdNucMap,
                                                        ipdGeneMap, descSeqMap)
if(debugging <= 3) { 
    err.println "${featureMap.size()} ids: " + featureMap.keySet().join(",")
}

// output
featureMap.each { k, vList ->
    //(id, gene, seqStartIndex, seqEndIndex) = parseDescription(desc)
    outFile = new File(k + ".ft.txt")
    FileWriter outWriter = new FileWriter(outFile)
    gbWriter = new PrintWriter(outFile.newOutputStream(), true)
    if(debugging <= 1) { 
        err.println "outputting to ${outFile}"
    }
    vList.each { v ->
        if(debugging <= 1) { 
            err.println "output feature ${k} ${v}"
        }
        outputGenBankSource(gbWriter, v)
        gbWriter.println ">Feature ${v.id}"
//todo        correctAugustus(v, ipdNucMap, ipdGeneMap) // reannotate 2DP1 and 3DP1
        //err.println "after correctAugustus full match, gene now ${v.gene}, pcall now ${v.pcall}, ccall now ${v.ccall}, gcall now ${v.gcall}"//todo
        outputGenBank(v, ipdNucMap, gbWriter)
    }
    gbWriter.close()
} // each feature
outputGLs(featureMap, outGLFile)
// end main

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
    e.PHASE = 6
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
 * @return a Map from ID to List of Expandos, each containing:
 *   id: id of the individual
 *   gene; the IPD name of the gene
 *   desc: the description of the sequence from the fasta
 *   sequence; the full sequence from the fasta as DNASequence
 *   startSeq: 5' start of the input sequence (1-based index)
 *   startSeq: 3' end of the input sequence (1-based index)
 *   start5p: 5' start of the KIR gene (1-based index)
 *   end3p: 3' end of the KIR gene (1-based index)
 *   pcall: the ultimate protein call/interpretation as a gl string (=~ 'NEW' if new)
 *   ccall: the ultimate cDNA call/interpretation as a gl string (=~ 'NEW' if new)
 *   gcall: the ultimate dna full gene call/interpretation as a gl string (=~ 'NEW' if new)
 *   exonStartTreeSet = [:] // exon index -> start of the exon (1-based index)
 *   exonEndTreeSet = [:] // feature -> end of the exon (1-based index)
 */
def Map<String, ArrayList<Expando>> processGFF(File gffFile, Map<String,String> ipdProtMap,
                                               Map<String,String> ipdNucMap,
                                               Map<String,String> ipdGeneMap,
                                               LinkedHashMap<String, DNASequence> descSeqMap) {
    if(debugging > 3) {
        err.println "processGFF(gffFile=${gffFile})"
    }

    Expando eFormat = initGFF() // properties are column properties

    // this will be returned; ID -> Expando
    Map <String, ArrayList<Expando>> geneMap = new HashMap()
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
        if(type != "transcription_start_site") {
            continue;
        }
        Expando eGene = processGFFGene(eFormat, line, listReader, descSeqMap, gffAAMap,
                                       gffCDSMap, gffCodingMap, ipdProtMap, ipdNucMap, ipdGeneMap)
        previousEntryList = geneMap[eGene.id]
        ArrayList<Expando> toAdd = new ArrayList()
        ArrayList<Expando> toRemove = new ArrayList()
        if(previousEntryList != null) { 
            previousEntryList.each { previousEntry ->
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
            ArrayList<Expando> al = new ArrayList()
            if(debugging <= 2) { 
                err.println "processGFF: adding ${eGene.id} " + eGene
            }
            al.add(eGene)
            geneMap[eGene.id] = al
        }
    } // each row

    if(debugging <= 3) { 
        err.println "processGFF: return ${geneMap.keySet().size()} IDs"
    }
    return geneMap
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
 * @return an Expando; see processGFF
 * @todo modularize this (see exon scope)
 */
def Expando processGFFGene(Expando eFormat, ArrayList line,
                           CsvListReader listReader, LinkedHashMap<String, DNASequence> descSeqMap,
                           LinkedHashMap<String, ProteinSequence> gffAAMap,
                           LinkedHashMap<String, DNASequence> gffCDSMap,
                           LinkedHashMap<String, DNASequence> gffCodingMap,
                           Map<String,String> ipdProtMap,
                           Map<String,String> ipdNucMap,
                           Map<String,String> ipdGeneMap) {
    e = new Expando() // return value
    e.exonStartTreeSet = new TreeSet()
    e.exonEndTreeSet = new TreeSet()
    HashSet<String> pCallSet = new HashSet()
    HashSet<String> cCallSet = new HashSet()
    HashSet<String> gCallSet = new HashSet()
    String gene // GFF gene name (e.g., 'g1')
    String transscript // GFF transcript name (e.g., 't1')
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
        if(type == "transcription_start_site") { // this precedes 'tss'
            // get the ID, and location of the sequence within the haplotype
            // e.g., cB02~tB01_GU182355.1_2DL4_67155-80842
            desc = line[eFormat.DESC]
            e.desc = desc
            (id, gene, seqStartIndex, seqEndIndex) = parseDescription(desc)
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
            e.startSeq = seqStartIndex
            e.endSeq = seqEndIndex
            
            attributes = line[eFormat.ATTRIBUTES]
            if(debugging <= 2) {
                err.println "processGFFGene: attributes=${attributes}"
            }
            // e.g., ID=g1.t1.tss1;Parent=g1.t1;
            (gene, transcript) = attributes.split('\\.')
            gene = gene.replaceFirst("ID=", "")
            ipdGene = GENE_NAME
            if(debugging <= 3) {
                err.println "processGFFGene: ${gene}.${transcript}"
                err.println "processGFFGene: ipdGene=${ipdGene}"
            }
            e.start5p = start
            e.end3p = end
            e.gene = GENE_NAME
        } else if(type == "CDS") { // or "exon"
            e.exonStartTreeSet.add(start)
            e.exonEndTreeSet.add(end)
            if(debugging <= 3) {
                err.println "processGFFGene: gene ${gene} exon ${exonIndex}: ${start}-${end}"
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
        } else if(type == "transcription_end_site") {
            if(debugging <= 2) {
                err.println "processGFFGene: found tts"
            }
            break; // time to returnfrom the method
        }
        // noop; see break at "tts"
        line = listReader.read()
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
        err.println "processGFFGene: return ${e.pcall}, ${e.ccall}, ${e.gcall}"
    }
    return e
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
    gffAAMap = FastaReaderHelper.readFastaProteinSequence(new File(protFileName))
    if(debugging <= 3) { 
        err.println gffAAMap.keySet().size() + " items"
    }

    // .cdsexons
    nucFileName = gffFile.getAbsolutePath().replace(".gff", ".cdsexons")
    if(debugging <= 3) { 
        err.println "loadGFFMaps: loading ${nucFileName}..."
    }
    gffCDSMap = FastaReaderHelper.readFastaDNASequence(new File(nucFileName))
    if(debugging <= 3) { 
        err.println gffCDSMap.keySet().size() + " items"
    }
    
    // .codingseq
    nucFileName = gffFile.getAbsolutePath().replace(".gff", ".codingseq")
    if(debugging <= 3) { 
        err.println "loadGFFMaps: loading ${nucFileName}..."
    }
    gffCodingMap = FastaReaderHelper.readFastaDNASequence(new File(nucFileName))
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
 * @param featureMap a Map from IPD/IMGT gene name to Expando (see processGFF())
 * @param outGLFile the file to which to write the GL strings
 * @todo rename: not gl any more
 */
def void outputGLs(Map<String, Expando> featureMap, File outGLFile) {
    if(debugging <= 1) {
        err.println "outputGLs(outGLFile=${outGLFile.getName()})"
    }
    FileWriter outGLWriter = new FileWriter(outGLFile)

    outGLWriter.println "ID\tgene\tprotein\tcDNA\tfull"
    featureMap.each { id, eList ->
        eList.each { e ->
            outGLWriter.println "${e.id}\t${e.gene}\t${e.pcall}\t${e.ccall}\t${e.gcall}"
        }
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
    outArg = new FileArgument("o", "ft", "output feature table file",
                              false)
    outGLArg = new FileArgument("l", "gl", "output gl file", true)
    geneSysArg = new StringArgument("g", "gene", "e.g., 'KIR2DL4' or 'HLA-A'",
                                    true)
    ArgumentList arguments = new ArgumentList(help, ipdArg, gffArg, seqArg, outArg,
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

    return [ipdArg.getValue(), gffArg.getValue(), seqArg.getValue(), outArg.getValue(),
            outGLArg.getValue(), geneSysArg.getValue()]
} // handleArgs

/*
 * outputGenBank
 * 
 * 1 based indexes
 * http://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html
 * this doesn't work for new protein level alleles (todo)
 * @param ipdNucMap Maps of description -> sequence for '_nuc.fasta' file from IPD-KIR
 */
def void outputGenBank(Expando query, Map<String,String> ipdNucMap, PrintWriter gbOut) {
    if(debugging <= 1) {
        err.println "outputGenBank(${query.gene})"
    }
    Integer startIndex = new Integer(query.startSeq)
    // gene
    gene = query.gene
    pseudo = false
    if(gene =~ /[23]DP/) {
        pseudo = true
    }
    geneStart = query.start5p // from 1-based to 0-based
    geneEnd = query.end3p
    if(startIndex != 1) {
        geneStart += startIndex
        geneEnd += startIndex
    }
    gbOut.println "${geneStart}\t${geneEnd}\tgene"
    gbOut.println "\t\t\tgene\t${query.gene}"
    outcall = query.gcall
    if(query.gcall =~ /${NEW_ALLELE}/) {
        outcall = query.ccall ? query.ccall : query.pcall
    }
    // todo: add "KIR" to call?
    gbOut.println "\t\t\tallele\t${outcall}"
    if(query.gcall =~ /${NEW_ALLELE}/) {
        gbOut.println "\t\t\tnote\tnew full-gene allele in ${NOMEN_VER}"
    }
    if(pseudo == true) {
        gbOut.println "\t\t\tnote\tpseudo"
        gbOut.println "\t\t\tnote\tpseudogene=\"unprocessed\""

    }
    gbOut.println "\t\t\tdb_xref\tinternal:${query.id}"
    
    // mRNA ('exon' in GFF)
    i = 0
    numExons = query.exonStartTreeSet.size()
    if(debugging <= 2) {
        err.println "outputGenBank: start mRNA"
        err.println "outputGenBank: ${numExons} exons"
    }
    if((numExons != 7) && (numExons != 8)) {
        err.println "WARNING ${numExons} exons for ${query.id}/${query.gene}"
    }
    // e.g., startSeq=67155, start5p=946, exonStartTreeSet=(0=981, ...), exonEndTreeSet=(0=1020, ...)
    // e.g., 5' UTR (start 5p - exon 1): 68101-68136  (68136-68101+1=36 bp)
    // e.g., exon 1 (index 0): 68136-68175  (68175-68136+1=40 bp)
    // e.g., exon 3 (index 2): 69277-69561   (69561-69277+1=285 bp)
    Iterator siter = query.exonStartTreeSet.iterator()
    Iterator eiter = query.exonEndTreeSet.iterator()
    for(i=0; siter.hasNext(); i++) { 
        //todo(remove) for(i=0; i < numExons; i++) {
        if(debugging <= 2) {
            err.println "outputGenBank: exon ${i}"
        }
        //start = query.exonStartMap[i] ? query.exonStartMap[i] : 0
        start = siter.next()
        if(i == 0) { // for exon 1, include 5' UTR
            if(query.start5p < start) { 
                start = query.start5p
            }
        }
        //todo end = query.exonEndMap[i] ? query.exonEndMap[i] : 0
        end = eiter.next()
/*        if((i+1) == numExons) { // for last exon, include 3' UTR
            end = query.end3p
        }
issue with e.g., 2DP1
*/
        length = end - start + 1
        if(startIndex != 1) {
            start += startIndex
            end += startIndex
        }
        if(debugging <= 2) {
            err.println "outputGenBank: ${start}-${end} = ${length}"
        }
pseudo
        if(pseudo == true) {
            gbOut.println "${start}\t${end}\texon"
            gbOut.println "\t\t\tgene=${query.gene}"
            gbOut.println "\t\t\tnumber=${i+1}"
        } else {
            if(i == 0) {
                gbOut.println "${start}\t${end}\tmRNA"
            } else { 
                gbOut.println "${start}\t${end}\t"
            }
        }
    } // each intron/exon pair
    gbOut.println "\t\t\tgene\t${query.gene}"
    if(GENE_SYSTEM == "KIR") { // todo: add this for HLA too
        outputProduct(query.gene, gbOut)
    }
    
    // CDS ('cds' in gff)
    if(!(query.gene =~ /\dDP/)) { // skip the pseudo genes
        // output the exon regions containing the coding sequence
        if(debugging <= 2) {
            err.println "outputGenBank: start CDS"
        }
        outputLength = 0
        Iterator siter2 = query.exonStartTreeSet.iterator()
        Iterator eiter2 = query.exonEndTreeSet.iterator()
        for(i=0; siter2.hasNext(); i++) { 
            start = siter2.next()
            if(start == null) {
                continue
            }
            end = eiter2.next()
            if(startIndex != 1) {
                start += startIndex
                end += startIndex
            }
            if(i == 0) {
                gbOut.println "${start}\t${end}\tCDS"
            } else { 
                gbOut.println "${start}\t${end}\t"
            }
            Integer length = end - start + 1
            // warn about KIR lengths here
            if(GENE_SYSTEM == "KIR") { // todo: implement for HLA too
                if(!KIR_EXON_SIZES.contains(length)) {
                    if(debugging <= 4) { 
                        err.println "WARNING: exon ${i} length ${length} (${start}-${end}) for ${query.id}/${query.gene}"
                    }
                }
            }
        } // each exon

        gbOut.println "\t\t\tgene\t${query.gene}"
        if(GENE_SYSTEM == "KIR") { // todo: implement for HLA too
            outputProduct(query.gene, gbOut)
        }
        gbOut.println "\t\t\tcodon_start\t1"
    } // CDS

    if(debugging <= 1) {
        err.println "outputGenBank: return"
    }
    return
} // outputGenBank

/*
 * Outputs the 'source' and 'misc_feature' GenBank sections.
 * @todo move this to combineTables
 */
def void outputGenBankSource(PrintWriter gbOut, Expando query) {
    if(debugging > 5) {
        err.println "outputGenBankSource()"
    }    
    Integer startIndex = 1
    Integer endIndex =  query.sequence.getLength() + 1
    // source
    gbOut.println "${startIndex}\t${endIndex}\tsource"
    gbOut.println "\t\t\torganism\tHomo sapiens"
    gbOut.println "\t\t\tmol_type=genomic DNA"
    gbOut.println "\t\t\tdb_xref\ttaxon:9606"
    chr = "6" // default to hla
    location = "6p21.3"
    if(GENE_SYSTEM == "KIR") {
        chr = "19"
        location = "19q13.4"
    }
    gbOut.println "\t\t\tchromosome\t${chr}"
    gbOut.println "\t\t\tmap\t${location}"

    // misc_feature
    gbOut.println "${startIndex}\t${endIndex}\tmisc_feature"
    gbOut.println "\t\t\tnote\t${GENE_SYSTEM} gene cluster"

    if(debugging > 5) {
        err.println "outputGenBankSource: return"
    }    
    return
} // outputGenBankSource

def void outputProduct(String gene, PrintWriter gbOut) {
    if(debugging <= 1) {
        err.println "outputProduct(gene=${gene})"
    }    

    gene = gene.replace("KIR", "").replace("HLA-", "")
    numDomains = gene[0].toInteger()
    tailLength = gene[2]
    number = gene[3].toInteger()

    domains = "two"
    length = "short cytoplasmic tail"
    if(numDomains == 3) {
        domains = "three"
    }
    if(tailLength == "L") {
        length = "long cytoplasmic tail"
    } else if(tailLength == "P") {
        length = "pseudo gene"
    }
    

    gbOut.println "\t\t\tproduct\tkiller cell immunoglobulin-like receptor ${domains} domains ${length} ${number}"
    if(debugging <= 1) {
        err.println "outputProduct: return"
    }    
    return
} // outputProduct

// todo: document
def Boolean sequenceEquals(org.biojava.nbio.core.sequence.template.Sequence a,
                           org.biojava.nbio.core.sequence.template.Sequence b) {
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
/*    if(!gene.contains("2DP1") && !gene.contains("3DP1")) {
        return
    }
*/
    cDNANew = new String() // the new cDNA sequence
    numExons = v.exonStartTreeSet.size()
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
