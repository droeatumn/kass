#!//usr/local/sdkman/candidates/groovy/current/bin/groovy
//todo #!/usr/bin/env groovy
/*
 * Orient sequences in a fasta.
 *
 * e.g., orient.groovy -i example1.contigs.fasta -p cap.fasta -o out.fasta
 *
 * Requires 
 *   - BioJava 4's Core jar: http://www.biojava/org
 *   e.g., 
 *     export CLASSPATH=$HOME/bin/jars/biojava4-core.jar:$CLASSPATH
 *
 *
 * @author Dave Roe
 */

import org.biojava.nbio.core.sequence.*
import org.biojava.nbio.core.sequence.io.*
import org.biojava.nbio.core.sequence.compound.*
import groovy.cli.commons.OptionAccessor
import groovy.cli.commons.CliBuilder

// things that may change per run
debugging = 3 // TRACE=1, WARN=2, DEBUG=3, INFO=4, ERROR=5

// things that probably won't change per run
err = System.err
probeOrder = [2, 3, 4, 7, 10, 12, 13, 1, 5, 6, 8, 9, 11, 14, 15, 16, 17, 18]

OptionAccessor options = handleArgs(args)

ambigDNA = new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet())

// descSeqMap: description -> DNASequence
FastaReader<DNASequence, NucleotideCompound> parentReader = new FastaReader<DNASequence, NucleotideCompound>(new File(options.i), new PlainFastaHeaderParser<DNASequence, NucleotideCompound>(), ambigDNA)
LinkedHashMap<String, DNASequence> descSeqMap = parentReader.process()
if(debugging <= 3) { 
    err.println "${descSeqMap.keySet().size()} input descriptions: " + descSeqMap.keySet().join(",")
}

// probeSeqMap: description -> DNASequence
parentReader = new FastaReader<DNASequence, NucleotideCompound>(new File(options.p), new PlainFastaHeaderParser<DNASequence, NucleotideCompound>(), ambigDNA)
LinkedHashMap<String, DNASequence> probeSeqMap = parentReader.process()
if(debugging <= 3) { 
    err.println "${probeSeqMap.keySet().size()} probe descriptions"
}

descSeqMap.each { desc, dnaSeq ->
    reverseIt = checkOrientation(dnaSeq, probeSeqMap)
    if(reverseIt) {
        dnaSeqNew = new DNASequence(dnaSeq.getReverseComplement().getSequenceAsString())
        dnaSeqNew.setOriginalHeader(dnaSeq.getOriginalHeader())
        descSeqMap[desc] = dnaSeqNew
    }
} // each sequence in the fasta

FastaWriterHelper.writeNucleotideSequence(new File(options.o), descSeqMap.values())

/*writer = new FastaWriterHelper(new File(options.o), descSeqMap.values(),
                               new GenericFastaHeaderFormat(), 1000000)
writer.process()*/

// end main

/*
 * checkOrientation
 */
Boolean checkOrientation(DNASequence seq,
                         LinkedHashMap<String, DNASequence> probeSeqMap) {
    Boolean ret = false
    for(s in probeSeqMap) {
        String desc = s.key
        DNASequence probeSeq = s.value

        probeStr = probeSeq.getSequenceAsString()
        seqStr = seq.getSequenceAsString()

        if(seqStr.contains(probeStr)) {
            break
        } else if(seqStr.contains(probeSeq.getReverseComplement().getSequenceAsString())) {
            ret = true
            break
        }
    } // each probe

    return ret
} // checkOrientation

/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'orient.groovy [options] ',
      header:'Options:')
    cli.help('print this message')
    cli.i(longOpt:'in', args:1, argName:'input', 'input FASTA file',
      required: true)
    cli.p(longOpt:'probe', args:1, argName:'probe', 'input probe file',
      required: true)
    cli.o(longOpt:'out', args:1, argName:'output', 'output FASTA file',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
