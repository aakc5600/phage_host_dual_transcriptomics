process GSE_TO_SRA {

    conda "${params.conda_path}/SRATOOLS"

    input :
    val geo

    output:
    path "**.sra", emit: sras

    script:
    """
    prefetch $geo
    """

    stub:
    """
    touch placeholder.sra placeholder2.sra
    """
}

process SRA_TO_FASTQ_SE {

    conda "${params.conda_path}/SRATOOLS"

    tag "fetching $params.geo"
    
    input:
    each path(sra)

    output:
    tuple env('SAMPLE'), path("*[0-9][0-9].sra.fastq.gz"), emit: singles
    
    script:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    fasterq-dump $sra
    find . -name "*.fastq" -exec sh -c 'pigz {} || gzip {}' \\;
    """

    stub:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}    
    touch \${SAMPLE}_STUB_15.fastq.gz
    """
}

process SRA_TO_FASTQ_PE {

    conda "${params.conda_path}/SRATOOLS"

    tag "fetching $params.geo"
    
    input:
    each path(sra)

    output:
    tuple env('SAMPLE'), path("*_{1,2}.sra.fastq.gz"), emit: pairs
    
    script:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    fasterq-dump $sra
    find . -name "*.fastq" -exec sh -c 'pigz {} || gzip {}' \\;
    """

    stub:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}    
    touch \${SAMPLE}_STUB_1.fastq.gz
    touch \${SAMPLE}_STUB_2.fastq.gz
    """
}


workflow FETCHREADS {
    take:
    geo
    reads
    pairedEnd

    main:
    if ( geo != "false" ) {
        sra = GSE_TO_SRA(geo)
        //singleSRAs = sra.sras.flatten() //flatMap {n -> [n[0], n[1]]}
        //singleSRAs.view()
        if ( pairedEnd == "true" ) {
            fastq = SRA_TO_FASTQ_PE(sra.sras)
            reads = fastq.pairs
            fastq.pairs.view()
        }
        else {
            fastq = SRA_TO_FASTQ_SE(sra.sras)
            reads = fastq.singles
            fastq.singles.view()
        } 
    }
    else {
        if ( pairedEnd == "true") {
            reads = channel.fromFilePairs(reads, size: 2, checkIfExists: true)
        }
        else {
            reads = channel.fromFilePairs(reads, size: -1, checkIfExists: true)
        }
    }

    emit:
    read_pairs_ch = reads
}