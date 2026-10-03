// Performs annotation of phage genomes with pharokka
process ANNOTATE_PHA {

    conda "${params.conda_path}/pharokka"

    tag "Annotating $fasta, using Pharokka on $threads threads"

    input:
    path fasta
    path db
    val threads

    output:
    path "pharokka/*.gbk", emit: out

    script:
    """
    pharokka run -d $db -t $threads -i '$fasta' -o pharokka
    """
}

// Performs annotation of phage genomes with phold
process ANNOTATE_PHO {

    conda "${params.conda_path}/pholdENV"

    tag "Annotating $gbk, using Phold on $threads threads"

    input:
    path gbk
    path db
    val threads

    output:
    path "phold/*.gbk", emit: out

    script:
    """
    phold run -d $db -t $threads --i '$gbk' -o phold
    """
}

// Performs annotation of phage genomes with phynteny
process ANNOTATE_PHY {

    conda "${params.conda_path}/phynteny_transformer"
    
    tag "Annotating $gbk, using Phynteny"

    input:
    path gbk
    path model

    output:
    path "phynteny/*.gbk", emit: out

    script:
    """
    phynteny_transformer -m $model -o phynteny '$gbk'
    """
}

// Performs Downstream processing automatically and includes annotation by phx tools
process AUTO_PROCESSING {
    conda "${params.conda_path}/DOWNSTREAM"

    tag "Processing $bulkPath, using $sampleDict"

    input:
    path bulkPath
    path metaPath
    path gffPath
    val sampleDict
    val middle
    val late
    path pwd

    output:
    path "*.tsv", emit: analysis

    script:
    """
    python "${pwd}/downstream_processing/downstream.py" --bulkPath $bulkPath --metaPath $metaPath --gffPath $gffPath --sampleDict $sampleDict --middle $middle --late $late
    """

    stub:
    """
    touch test.tsv
    """
}

process AUTO_PROCESSING_PH {
    conda "${params.conda_path}/DOWNSTREAM"

    tag "Processing $bulkPath, using $sampleDict and $gbkPath"

    input:
    path bulkPath
    path metaPath
    path gffPath
    path gbkPath
    val sampleDict
    val middle
    val late
    path pwd

    output:
    path "*.tsv", emit: analysis

    script:
    """
    python "${pwd}/downstream_processing/downstream.py" --bulkPath $bulkPath --metaPath $metaPath --gffPath $gffPath --gbkPath $gbkPath --sampleDict $sampleDict --middle $middle --late $late
    """

    stub:
    """
    touch test.tsv
    """
}
