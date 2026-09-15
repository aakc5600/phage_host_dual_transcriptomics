process TRIM_SE {

    conda "${params.conda_path}/RNASEQ"

    tag "Trimming single-end reads via cutadapt $sampleID"
    //publishDir "$params.outputDir/cutadapt", pattern: "*.cutadapt.json", mode: params.pubDirMode

    input:
    tuple val(sampleID), path(reads)

    output:
    val sampleID, emit: sampleID
    path "${sampleID}_trimmed.fastq", emit: trimmedReads
    path "${sampleID}.cutadapt.json", emit: cutadaptReport

    script:
    """
    cutadapt -a $params.adapter1 -g $params.adapter2 --json=${sampleID}.cutadapt.json -q 28 --quality-base 33 -o ${sampleID}_trimmed.fastq $reads
    """
}

process TRIM_PE {

    conda "${params.conda_path}/RNASEQ"

    tag "Trimming paired-end reads via cutadapt $sampleID"
    //publishDir "$params.outputDir/cutadapt", pattern: "*.cutadapt.json", mode: params.pubDirMode

    input:
    tuple val(sampleID), path(reads)

    output:
    val sampleID, emit: sampleID
    path "${sampleID}_trimmed*.fastq", emit: trimmedReads
    path "${sampleID}.cutadapt.json", emit: cutadaptReport

    script:
    """
    cutadapt -a $params.adapter1 -A $params.adapter2 --json=${sampleID}.cutadapt.json -q 28 --quality-base 33 -o ${sampleID}_trimmed_R1.fastq -p ${sampleID}_trimmed_R2.fastq $reads
    """
}
