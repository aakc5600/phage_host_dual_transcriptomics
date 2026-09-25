def parseHelp(help, verbose) {
    if (verbose) {
        helpMessageVerbose()
        exit 0
    }
    if (help) {
        helpMessage()
        exit 0
    }
}

def helpMessage() {
    log.info"""
    ==-------------------------------------------------------------------==
    Pipeline for processing dual RNA-seq data for the PhageExpression Atlas
    ==-------------------------------------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run rnaseq_workflow.nf --reads '*_{1,2}.fastq.gz' --pairedEnd --hostGenome 'hostGenome.fasta' --hostGFF 'hostGenome.gff' --phageGenome 'phageGenome.fasta' --phageGFF 'phageGenome.gff3'

    Arguments:
      --reads                 Path specifying the reads. E.g. ./*_reads.fastq.gz for single-end, or ./*_reads_{1,2}.fastq.gz
      --geo                   Optional. Overrides --reads. Path specifying an accession number to fetch from GEO.

      --outputDir             Output directory.
      --hostGenome            Host genome path in fasta format.
      --phageGenome           Phage genome path in fasta format.
      --hostGFF               Host genome annotation path in gff3 format.
      --phageGFF              Phage genome annotation path in gff3 format.

      --pairedEnd             Default is false. Set, if paired-end reads are supplied.

      --adapter1              Default is sequence of TruSeq adapter.
      --adapter2              Default is sequence of TruSeq adapter.
    
      --conda_path            Default is "/home/$params.user/miniconda3/envs"

    Arguments can also be set in the config file (./conf/params.config)
    """.stripIndent()
}

def helpMessageVerbose() {
    log.info"""
    ==-------------------------------------------------------------------==
    Pipeline for processing dual RNA-seq data for the PhageExpression Atlas
    ==-------------------------------------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run rnaseq_workflow.nf --reads '*_{1,2}.fastq.gz' --pairedEnd --hostGenome 'hostGenome.fasta' --hostGFF 'hostGenome.gff' --phageGenome 'phageGenome.fasta' --phageGFF 'phageGenome.gff3'

    Arguments:
      --reads                 Path specifying the reads. E.g. ./*_reads.fastq.gz for single-end, or ./*_reads_{1,2}.fastq.gz Default is "/home/$params.user/input/reads/*.fast*"
        [symbols in {} will be used to match paired end read files together. Paired filenames for paired reads must be identical before the {}]
        [Both gzip compressed and raw fastq files can be used]
        
      --geo                   Optional. Overrides --reads. Path specifying an accession number to fetch from GEO.
      [Supports both single-end and paired-end data]
      [Requires ncbi sra-tools to be installed]

      --outputDir             Output directory.
      --inputDir              Optional. Sets default path for reads and genomes. Default is "/home/$params.user/input"
        [only used if no path for reads or genomes/gffs is provided]
      --hostGenome            Host genome path in fasta format. Default/Currently set as "/home/$params.user/input/genomes/host.fna"
      --phageGenome           Phage genome path in fasta format. Default/Currently set as "/home/$params.user/input/genomes/phage.fna"
      --hostGFF               Host genome annotation path in gff3 format. Default/Currently set as "/home/$params.user/input/genomes/host.gff"
      --phageGFF              Phage genome annotation path in gff3 format. Default/Currently set as "/home/$params.user/input/genomes/phage.gff"

      --pairedEnd             Default is false. Set, if paired-end reads are supplied.
        [If set, provided filenames for reads must contain {X,Y} (see --reads above)]

      --adapter1              Default is sequence of TruSeq adapter.
        [Default AGATCGGAAGAGCACACGTCTGAACTCCAGTCA]
      --adapter2              Default is sequence of TruSeq adapter.
        [Default AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]

      --conda_path            Default is "/home/$params.user/miniconda3/envs"
        [Path to environments folder for your install of conda]

      DOWNSTREAM PROCESSING
      --auto                  If true, downstream analysis will be automatically started
      [requires metaPath, sampleDict and either timeStart and timeEnd or timeUntilMiddle and timeUntilLate to be set]
      [Samples are treated as single dataset. Different treatments (eg. +- antibiotic) are not automatically parsed]
      [If further parsing of data is required, downstream analysis has to be carried out manually via the provided jupyter notebook]

      --sampleDict            Dictionary-like string denoting which SRR number corresponds to which sample for downstream analysis
      [Format must be '[SRR-number]:[timepoint_reapeat], [SRR-number]:[timepoint_reapeat], etc.'     eg. see ./conf/params.config]

      --metaPath              sraRunTable.csv to use for downstream processing. Default/Currently set as "/home/$params.user/input/metadata/sraRunTable.csv"

      --timeStart             Timepoint [min] post infection of first sample.
      --timeEnd               Timepoint [min] post infection of last sample.
      [Will only be used to provide default values if timeUntilMiddle/Late is not explicitly set]
      --timeUntilMiddle       Timepoints [min] after which genes are counted as middle genes (based on maximal or threshold of 20% expression).
      --timeUntilLate         Timepoints [min] after which genes are counted as late genes (based on maximal or threshold of 20% expression)
      [It is reccommended to set timeUntilXXXXs explicitly rather than relying on default values inferred by timeStart and timeEnd (1/4 and 1/2 of the difference between those)]
      Arguments can also be set in the config file (./conf/params.config)
    """.stripIndent()
}