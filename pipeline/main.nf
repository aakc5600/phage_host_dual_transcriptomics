#!/usr/bin/env nextflow

include { parseHelp; helpMessage; helpMessageVerbose } from './modules/helpMsg'
include { FETCHREADS } from './subworkflows/fetchReads'
include { BUILDHISAT2BASE } from './modules/buildHisatBase'
include { PROCESSRNASEQ } from './subworkflows/perReadProcessing'
include { FEATURECOUNTS } from './modules/readCounting'
include { CLEANUP } from './modules/cleanDir'
include { MULTIQC } from './modules/multiQC'

workflow {
    main:
    parseHelp(params.help, params.verboseHelp)    
    hisat2 = BUILDHISAT2BASE(params.hostGenome, params.phageGenome, params.hostGFF, params.phageGFF)
    read_pairs_ch = FETCHREADS(params.geo, params.reads, params.pigz)
    rnaseq = PROCESSRNASEQ(read_pairs_ch, BUILDHISAT2BASE.out.alignmentBase)
    featurecounts = FEATURECOUNTS(PROCESSRNASEQ.out.sortedBamFile.collect(), BUILDHISAT2BASE.out.dualGFF)
    multiqc = MULTIQC(FEATURECOUNTS.out.countSummary, PROCESSRNASEQ.out.cutadaptReport.collect().ifEmpty([]), PROCESSRNASEQ.out.fastqc_pre.collect(), PROCESSRNASEQ.out.fastqc_post.collect())    

    publish:
    // This looks so ugly, but idk if there's a nicer way 
    reads = read_pairs_ch.read_pairs_ch

    fastqc_pre = rnaseq.fastqc_pre
    fastqc_post = rnaseq.fastqc_post
    sortedBamFile = rnaseq.sortedBamFile
    bamIndex = rnaseq.bamIndex
    cutadaptReport = rnaseq.cutadaptReport

    alignmentBase = hisat2.alignmentBase
    dualGFF = hisat2.dualGFF

    countTable = featurecounts.countTable
    countSummary = featurecounts.countSummary

    multiqc_r = multiqc.multiqc_r
}
// pooling the outputs by module/subworkflow of origin
output {
    reads{ path "reads" }
    fastqc_pre{ path "rnaseq" }
    fastqc_post{ path "rnaseq" }
    sortedBamFile{ path "rnaseq" }
    bamIndex{ path "rnaseq" }
    cutadaptReport{ path "rnaseq" }

    alignmentBase{ path "hisat2" }
    dualGFF{ path "hisat2" }

    countTable{ path "featurecounts"}
    countSummary{ path "featurecounts"}

    multiqc_r{ path "multiqc"}
}