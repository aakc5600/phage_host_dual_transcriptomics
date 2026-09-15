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
      --hostGenome            Host genome path in fasta format.
      --phageGenome           Phage genome path in fasta format.
      --hostGFF               Host genome annotation path in gff3 format.
      --phageGFF              Phage genome annotation path in gff3 format.

      --pairedEnd             Default is false. Set, if paired-end reads are supplied.

      --adapter1              Default is sequence of TruSeq adapter.
      --adapter2              Default is sequence of TruSeq adapter.
    
      --conda_path            Default is "/home/$params.user/miniconda3/envs"

    Arguments can also be set in the config file (./conf/params.config)
    """.stripIndent()
}

def helpMessageVerbose() {
    log.info"""
    ==-------------------------------------------------------------------==
    Pipeline for processing dual RNA-seq data for the PhageExpression Atlas
    ==-------------------------------------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run rnaseq_workflow.nf --reads '*_{1,2}.fastq.gz' --pairedEnd --hostGenome 'hostGenome.fasta' --hostGFF 'hostGenome.gff' --phageGenome 'phageGenome.fasta' --phageGFF 'phageGenome.gff3'

    Arguments:
      --reads                 Path specifying the reads. E.g. ./*_reads.fastq.gz for single-end, or ./*_reads_{1,2}.fastq.gz Default is ""/home/$params.user/input/reads/*.fast*"
        [symbols in {} will be used to match paired end read files together. Paired filenames for paired reads must be identical before the {}]
        [Both gzip compressed and raw fastq files can be used]
      --outputDir             Output directory.
      --inputDir              Optional. Sets default path for reads and genomes. Default is "/home/$params.user/input"
        [only used if no path for reads or genomes/gffs is provided]
      --hostGenome            Host genome path in fasta format. Default is "/home/$params.user/input/genomes/host.fna"
      --phageGenome           Phage genome path in fasta format. Default is "/home/$params.user/input/genomes/phage.fna"
      --hostGFF               Host genome annotation path in gff3 format. Default is "/home/$params.user/input/genomes/host.gff"
      --phageGFF              Phage genome annotation path in gff3 format. Default is "/home/$params.user/input/genomes/phage.gff"

      --pairedEnd             Default is false. Set, if paired-end reads are supplied.
        [If set, you must also provide filenames must contain {X,Y} (see --reads above)]

      --adapter1              Default is sequence of TruSeq adapter.
        [Default AGATCGGAAGAGCACACGTCTGAACTCCAGTCA]
      --adapter2              Default is sequence of TruSeq adapter.
        [Default AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]

      --conda_path            Default is "/home/$params.user/miniconda3/envs"
        [Path to environments folder for your install of conda]

    Arguments can also be set in the config file (./conf/params.config)
    """.stripIndent()
}

workflow {
    main:
    if (params.help){
        if (params.verbose){
            helpMessageVerbose()
        }
        else {
            helpMessage()
        }
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