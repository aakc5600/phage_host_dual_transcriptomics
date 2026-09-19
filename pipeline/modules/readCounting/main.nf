

process FEATURECOUNTS {

    conda "${params.conda_path}/RNASEQ"

    tag "Count reads for all samples"
    //publishDir "$params.outputDir/countData", pattern: "*.tsv", mode: params.pubDirMode
    //publishDir "$params.outputDir/countData", pattern: "*.tsv.summary", mode: params.pubDirMode
    
    input:
    path sortedBamFiles
    path inputGFF

    output:
    path "countData.tsv", emit: countTable
    path "countData.tsv.summary", emit: countSummary    

    script:
    """
    featureCounts -p -M -O --primary \\
    -a $inputGFF \\
    -t $params.countFeature \\
    -o countData.tsv \\
    -g $params.featureIdentifier \\
    $sortedBamFiles
    """

    stub:
    """
    touch countData.tsv
    touch countData.tsv.summary
    """
}
