#!/usr/bin/env nextflow


include { BUILDHISAT2BASE } from './modules/buildHisatBase'
include { PROCESSRNASEQ } from './subworkflows/perReadProcessing'
include { FEATURECOUNTS } from './modules/readCounting'
include { CLEANUP } from './modules/cleanDir'
include { MULTIQC } from './modules/multiQC'


def helpMessage() {
    log.info"""
    ==-------------------------------------------------------------------==
    Pipeline for processing dual RNA-seq data for the PhageExpression Atlas
    ==-------------------------------------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run rnaseq_workflow.nf --reads '*_{1,2}.fastq.gz' --pairedEnd --hostGenome 'hostGenome.fasta' --hostGFF 'hostGenome.gff' --phageGenome 'phageGenome.fasta' --phageGFF 'phageGenome.gff3'

    Arguments:
      --reads                 Path specifying the reads. E.g. ./*_reads.fastq.gz for single-end, or ./*_reads_{1,2}.fastq.gz
      --outputDir             Output directory.
      --inputDir              Input directory, optional.
      --hostGenome            Host genome path in fasta format.
      --phageGenome           Phage genome path in fasta format.
      --hostGFF               Host genome annotation path in gff3 format.
      --phageGFF              Phage genome annotation path in gff3 format.

      --pairedEnd             Default is false. Set, if paired-end reads are supplied.

      --adapter1              Default is sequence of TruSeq adapter.
      --adapter2              Default is sequence of TruSeq adapter.

    Arguments can also be set in the config file (./conf/params.config)
    """.stripIndent()
}


workflow {
    main:
    if (params.help){
        helpMessage()
        exit 0
    }
    
    hisat2 = BUILDHISAT2BASE(params.hostGenome, params.phageGenome, params.hostGFF, params.phageGFF)
    rnaseq = PROCESSRNASEQ(params.reads, BUILDHISAT2BASE.out.alignmentBase)    
    featurecounts = FEATURECOUNTS(PROCESSRNASEQ.out.sortedBamFile.collect(), BUILDHISAT2BASE.out.dualGFF)
    multiqc = MULTIQC(FEATURECOUNTS.out.countSummary, PROCESSRNASEQ.out.cutadaptReport.collect().ifEmpty([]), PROCESSRNASEQ.out.fastqc_pre.collect(), PROCESSRNASEQ.out.fastqc_post.collect())    
    // grabbing all the module outputs
    // getting PROCESSRNASEQ() outputs
    fastqc_pre = rnaseq.fastqc_pre
    fastqc_post = rnaseq.fastqc_post
    sortedBamFile = rnaseq.sortedBamFile
    bamIndex = rnaseq.bamIndex
    cutadaptReport = rnaseq.cutadaptReport
    // getting BUILDHISAT2BASE() outputs
    alignmentBase = hisat2.alignmentBase
    dualGFF = hisat2.dualGFF
    // getting FEATURECOUNTS() outputs
    countTable = featurecounts.countTable
    countSummary = featurecounts.countSummary
    // getting MULTIQC() outputs
    multiqc_r = multiqc.multiqc_r

    publish:
    fastqc_pre = fastqc_pre
    fastqc_post = fastqc_post
    sortedBamFile = sortedBamFile
    bamIndex = bamIndex
    cutadaptReport = cutadaptReport

    alignmentBase = alignmentBase
    dualGFF = dualGFF

    countTable = countTable
    countSummary = countSummary

    multiqc_r = multiqc_r
}

output {
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