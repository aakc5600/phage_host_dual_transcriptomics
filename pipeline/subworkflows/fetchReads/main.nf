process fetch_reads_SE {
    input:
    each sample

    output:
    tuple val("${sample}"), path("${sample}.fastq.gz"), emit: read_list

    script:
    """
    fastq-dump --gzip "${sample}"
    """

    stub:
    """
    touch "${sample}.fastq.gz"
    """
}

process fetch_reads_PE {
    input:
    each sample

    output:
    tuple val("${sample}"), path("${sample}_{1,2}.fastq.gz"), emit: read_list

    script:
    """
    fastq-dump --gzip --split-3 "${sample}" 
    """

    stub:
    """
    touch "${sample}_1.fastq.gz" && touch "${sample}_2.fastq.gz"
    """    
}

workflow FETCHREADS {
    take:
    srrList
    reads

    main:




    if (params.srrList){
        samples = []
        count = 0
        file(srrList).eachLine { line ->
            samples[count] = line
            count += 1
        }
        if (params.pairedEnd) {
            reads = fetch_reads_PE(samples)
        }
        else {
            reads = fetch_reads_SE(samples)
        }
    }
    else {
        if (params.pairedEnd) {
            reads = channel.fromFilePairs(reads, size: 2, checkIfExists: true)
        }
        else {
            read_pairs_ch = channel.fromFilePairs(reads, size: -1, checkIfExists: true)
        }       
    }



    emit:
    read_pairs_ch = reads
}