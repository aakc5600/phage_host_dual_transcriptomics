/*
Filtering/processing alignment via samtools
*/

process FILTERSAMTOBAM {

    conda "${params.conda_path}/RNASEQ"

    tag "Processing alignment file ${sampleID}"
    //publishDir "$params.outputDir/alignments", pattern: "*.bam", mode: params.pubDirMode
    //publishDir "$params.outputDir/alignments", pattern: "*.bai", mode: params.pubDirMode

    input:
    val sampleID
    path samFile

    output:
    path "${sampleID}_sorted.bam", emit: sortedBamFile
    path "${sampleID}_sorted.bam.bai", emit: baiFile
    val sampleID, emit: sampleID

    script:
    """
    samtools view -bS $samFile > ${sampleID}.bam
    samtools sort ${sampleID}.bam -o ${sampleID}_sorted.bam
    samtools index ${sampleID}_sorted.bam  
    """

}


process BAMFORCOVERAGE {

    conda "${params.conda_path}/RNASEQ"

    tag "Processing bam file for coverage (wiggle files) for ${sampleID}"
    //publishDir "${params.outputDir}/alignments", mode: pubDirMode

    input:
    val sampleID
    path sortedBamFile

    output:
    path "${sampleID}_Chr_primary.bam", emit: filteredBamFile
    val sampleID, emit: sampleID

    script:
    """
    samtools view -b $sortedBamFile Chr > ${sampleID}_Chr.bam
    samtools view -F 256 ${sampleID}_Chr.bam > ${sampleID}_Chr_primary.bam
    """

}