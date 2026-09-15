process MULTIQC {

    conda "${params.conda_path}/MULTIQC"

    tag "multiQC analysis of processes"
    //publishDir "$params.outputDir/multiQC", pattern: "*.html", mode: params.pubDirMode

    input:
    path countSummary
    path cutadaptReport
    path fastqc_pre
    path fastqc_post

    output:
    path "multiQC_report.html", emit: multiqc_r 

    script:
    """
    multiqc . -n multiQC_report.html --no-data-dir -m custom_content -m preseq -m rseqc -m featurecounts -m star -m cutadapt -m fastqc -m qualimap -m salmon
    """
}
