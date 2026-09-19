/*process GSE_TO_SRA {
    input :
    val geo

    output:
    path "**.sra"

    script:
    """
    prefetch $geo
    """

    stub:
    """
    touch placeholder.sra
    """
}

process SRA_TO_FASTQ_SE {
    input:
    path sra

    output:
    tuple env('SAMPLE'), path("*.fastq.gz")

    script:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    fasterq-dump $sra
    gzip \$SAMPLE.fastq
    """

    stub:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    touch placeholder.fastq.gz
    """
}

process SRA_TO_FASTQ_PE {
    input:
    path sra

    output:
    tuple env('SAMPLE'), path("*_{1,2}.fastq.gz")

    script:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    fasterq-dump $sra
    find . -name '*.fastq' -exec gzip {} \\;
    """

    stub:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}    
    touch placeholder_1.fastq.gz
    touch placeholder_2.fastq.gz
    """
}

process SRA_TO_FASTQ_SE_PIGZ {
    input:
    path sra

    output:
    tuple env('SAMPLE'), path("*.fastq.gz")

    script:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    fasterq-dump $sra
    pigz \$SAMPLE.fastq
    """

    stub:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    touch placeholder.fastq.gz
    """
}

process SRA_TO_FASTQ_PE_PIGZ {
    input:
    path sra

    output:
    tuple env('SAMPLE'), path("*_{1,2}.fastq.gz")

    script:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}
    fasterq-dump $sra
    find . -name "*.fastq" -exec pigz {} \\;
    """

    stub:
    """
    SRA=$sra
    SAMPLE=\${SRA%%.sra}    
    touch placeholder_1.fastq.gz
    touch placeholder_2.fastq.gz
    """
}
*/

process SRA_TO_FASTQ {
    tag "fetching $geo"
    
    input:
    path sra

    output:
    tuple env('SAMPLE'), path("*[0-9][0-9].fastq.gz"), emit: singles, optional: true
    tuple env('SAMPLE'), path("*_{1,2}.fastq.gz"), emit: pairs, optional: true
    
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
    touch placeholder_1.fastq.gz
    touch placeholder_2.fastq.gz
    """
}


workflow FETCHREADS {
    take:
    geo
    reads
    // pigz

    main:
    if (params.geo) {
        sra = GSE_TO_SRA(params.geo)
        fastq = SRA_TO_FASTQ(sra)
        fastq.singles.view()
        fastq.pairs.view()
        if (params.pairedEnd) {
            reads = fastq.pairs
        }
        else {
            reads = fastq.singles
        }
        
        // if (params.pairedEnd) {
        //     if (pigz) {
        //         reads = SRA_TO_FASTQ_PE_PIGZ(sra)
        //     }
        //     else {
        //         reads = SRA_TO_FASTQ_PE(sra)
        //     }
        // }
        // else {
        //     if (pigz) {
        //         reads = SRA_TO_FASTQ_SE_PIGZ(sra)
        //     }
        //     else {
        //     reads = SRA_TO_FASTQ_SE(sra)
        //     }
        // }        
    }
    else {
        if (params.pairedEnd) {
            reads = channel.fromFilePairs(params.reads, size: 2, checkIfExists: true)
        }
        else {
            reads = channel.fromFilePairs(params.reads, size: -1, checkIfExists: true)
        }
    }

    emit:
    read_pairs_ch = reads
}