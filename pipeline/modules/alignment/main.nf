/*
    Performs alignment with hisat2
*/
process MAPPINGPE {

    conda "${params.conda_path}/RNASEQ"

    tag "Alignment of sample ${sampleID}"
    //publishDir "$params.outputDir/alignments", pattern: "*.hisat2.summary.log", mode: params.pubDirMode
    
    input:
    val sampleID
    path trimmedReads
    path alignmentBase

    output:
    val sampleID, emit: sampleID
    path "${sampleID}.sam", emit: samFile
    path "${sampleID}.hisat2.summary.log", emit: summaryFile

    script:
    """
    hisat2 -q -p 16 -x ${alignmentBase}/Hisat2Base \\
    -1 ${trimmedReads[0]} \\
    -2 ${trimmedReads[1]} \\
    --summary-file ${sampleID}.hisat2.summary.log \\
    -S ${sampleID}.sam
    """
}

process MAPPINGSE {

    conda "${params.conda_path}/RNASEQ"

    tag "Alignment of sample ${sampleID}"
    
    input:
    val sampleID
    path trimmedReads
    path alignmentBase

    output:
    val sampleID, emit: sampleID
    path "${sampleID}.sam", emit: samFile
    path "${sampleID}.hisat2.summary.log", emit: summaryFile

    script:
    """
    hisat2 -q -p 16 -x ${alignmentBase}/Hisat2Base \\
    -U $trimmedReads \\
    --summary-file ${sampleID}.hisat2.summary.log \\
    -S ${sampleID}.sam
    """
}