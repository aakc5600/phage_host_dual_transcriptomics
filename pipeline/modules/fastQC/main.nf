/*
    Performs fastqc analysis of a fastq file (gz or unzipped; untrimmed)
*/
process FASTQC_PE {

    conda "${params.conda_path}/RNASEQ"

    tag "FastQC analysis of sample ${sample_id}"
    //publishDir "$params.outputDir/fastQC", mode: params.pubDirMode
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_*", emit: fastqc_pre

    script:
    """
    mkdir fastqc_${sample_id}_1
    mkdir fastqc_${sample_id}_2
    fastqc -o fastqc_${sample_id}_1 -q ${reads[0]}
    fastqc -o fastqc_${sample_id}_2 -q ${reads[1]}
    """
}

process FASTQC_SE {

    conda "${params.conda_path}/RNASEQ"

    tag "FastQC analysis of sample ${sample_id}"
    //publishDir "$params.outputDir/fastQC", mode: params.pubDirMode
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_*", emit: fastqc_pre

    script:
    """
    mkdir fastqc_${sample_id}_1
    fastqc -o fastqc_${sample_id}_1 -q ${reads[0]}
    """
}