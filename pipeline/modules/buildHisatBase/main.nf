/* 
Building the HISAT2 base for alignment
*/

process BUILDHISAT2BASE {

    tag "Build dual reference genome and Hisat2Base for alignment"
    // publishDir "$params.outputDir/alignments", pattern: "*.ht2", mode: params.pubDirMode
    // publishDir "$params.outputDir/alignments", pattern: "*.gff3", mode: params.pubDirMode

    conda "${params.conda_path}/RNASEQ"

    input:
    path hostGenome
    path phageGenome
    path hostGFF
    path phageGFF

    output:
    path "Hisat2Base", emit: alignmentBase
    path "dualGenome.gff3", emit: dualGFF

    script:
    """
    cat $hostGenome $phageGenome > dualGenome.fasta
    cat $hostGFF $phageGFF > dualGenome.gff3
    mkdir Hisat2Base
    hisat2-build -f -p 16 dualGenome.fasta Hisat2Base/Hisat2Base
    """
}