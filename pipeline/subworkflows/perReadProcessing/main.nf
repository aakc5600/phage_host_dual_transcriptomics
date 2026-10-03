
include { FASTQC_PE } from '../../modules/fastQC'
include { FASTQC_SE } from '../../modules/fastQC'
include { FASTQCTRIMMED_PE } from '../../modules/fastQCTrimmed'
include { FASTQCTRIMMED_SE } from '../../modules/fastQCTrimmed'
include { MAPPINGPE } from '../../modules/alignment'
include { MAPPINGSE } from '../../modules/alignment'
include { FILTERSAMTOBAM } from '../../modules/filteringAlignment'
include { TRIM_SE } from '../../modules/trimming'
include { TRIM_PE } from '../../modules/trimming'

workflow PROCESSRNASEQ {
    take: 
        read_pairs_ch
        alignmentBase
        threads
        mem

    main: 
        if(params.pairedEnd) {
            // read_pairs_ch = channel.fromFilePairs(reads, size: 2, checkIfExists: true)
            fastqc_pre = FASTQC_PE(read_pairs_ch, threads)
            TRIM_PE(read_pairs_ch)
            fastqc_post = FASTQCTRIMMED_PE(TRIM_PE.out.sampleID, TRIM_PE.out.trimmedReads, threads)
            MAPPINGPE(TRIM_PE.out.sampleID, TRIM_PE.out.trimmedReads, alignmentBase)
            FILTERSAMTOBAM(MAPPINGPE.out.sampleID, MAPPINGPE.out.samFile, threads, mem)
            cutadaptReport = TRIM_PE.out.cutadaptReport
        }
        else {
            // read_pairs_ch = channel.fromFilePairs(reads, size: -1, checkIfExists: true)
            fastqc_pre = FASTQC_SE(read_pairs_ch, threads)
            TRIM_SE(read_pairs_ch)
            fastqc_post = FASTQCTRIMMED_SE(TRIM_SE.out.sampleID, TRIM_SE.out.trimmedReads, threads)
            MAPPINGSE(TRIM_SE.out.sampleID, TRIM_SE.out.trimmedReads, alignmentBase)
            FILTERSAMTOBAM(MAPPINGSE.out.sampleID, MAPPINGSE.out.samFile, threads, mem)
            cutadaptReport = TRIM_SE.out.cutadaptReport
        }

    emit:
        sortedBamFile = FILTERSAMTOBAM.out.sortedBamFile
        bamIndex = FILTERSAMTOBAM.out.baiFile
        cutadaptReport = cutadaptReport
        fastqc_pre = fastqc_pre
        fastqc_post = fastqc_post
}