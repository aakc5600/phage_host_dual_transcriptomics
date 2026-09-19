process AUTO_PROCESSING {

    conda "${params.conda_path}/RNASEQ"

    tag "Processing $bulkPath, using $sampleDict"

    input:
    path bulkPath
    path metaPath
    path gffPath
    val sampleDict
    path tools

    output:
    path "*.tsv", emit: analysis

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA
    from $tools import *
    import re


    # ## 1 Load the dataset and annotation
    # 
    # We start by loading the raw counts table, metadata and concatenated gff3 file including host and phage.

    # Raw counts table, output from nextflow pipeline
    bulkPath = '$bulkPath' 
    # Metadata from SRA. Can be aquired from SRA Run Selector
    metaPath = '$metaPath'
    # Concatenated phage/host gff3 file, output from nextflow pipeline
    gffPath = '$gffPath'

    # Load data as pandas dataframe
    df_initial = pd.read_csv(bulkPath, sep = '\\t', comment='#', index_col=0)
    metadata = pd.read_csv(metaPath)

    # ## 2 Format the dataset
    # ### 2.1 Annotate sample names# 
    # Issue now is that the metadata does not properly annoate samples names. We do this manually by adding another sample name column to the metadata with proper sample annotation (timepoint + replicate) based on GEO Sample accessions (GSExyz) using `sampleDict`.
    # # `annotateData(metadata, sampleDict)` adds the new sample name column and indexes the dataframe with SRA Run accession number matching the counts table column names (SRRxyz_sorted.bam) to prepare updating counts table column names.

    # Match GSM IDs and SampleNames inferred from GEO
    sampleDict = {'GSM6447614': '0_R1', 'GSM6447615': '0_R2', 'GSM6447616': '0_R3', 'GSM6447617': '1_R1', 'GSM6447618': '1_R2', 'GSM6447619': '1_R3',
                'GSM6447620': '4_R1', 'GSM6447621': '4_R2', 'GSM6447622': '4_R3', 'GSM6447623': '7_R1', 'GSM6447624': '7_R2', 'GSM6447625': '7_R3',
                'GSM6447626': '20_R1', 'GSM6447627': '20_R2', 'GSM6447628': '20_R3'}

    sampleDict = parseDict(${sampleDict})

    # Update metadata dataframe
    metadataFull = annotateData(metadata, sampleDict)

    # Based on the updated metadata, we can now rename the columns of our counts table by running `changeColnames`. Input is:
    # 
    # 1. Raw counts table without gene Chr, Start, End, Strand, Length information `df_initial.iloc[:,5:df_initial.shape[1]]`
    # 2. and updated metadata `metadataFull`
    # 
    # Finally, we can sort the columns properly.

    # Update column names of raw counts table.
    df = changeColnames(df_initial.iloc[:,5:df_initial.shape[1]], metadataFull)

    # Sort columns properly
    # allows processing of partial data for testing and simplifies usage. However, requires column names to have similar format to 'timepoint_replicate'
    df = df.sort_index(axis=1)
    # df = df[['0_R1', '0_R2', '0_R3', '1_R1', '1_R2', '1_R3', '4_R1', '4_R2', '4_R3', '7_R1', '7_R2', '7_R3', '20_R1', '20_R2', '20_R3']]

    # ### 2.2 Remove rRNA genes
    # 
    # To perform depletion of ribosomal RNA genes, which quantitatively make up the majority of RNA within a sample, we consult the gff3 file with host and phage information to get rRNA GeneIDs. We do this by filtering the gff3 file for gene entries and accessing the attributes column to get GeneIDs, GeneType and GeneSymbol (important for later processing). Example:
    # 
    # - GeneID: gene-b0001
    # - GeneType: protein_coding
    # - GeneSymbol: thrL (Important note: gene symbols not available for most genes.)
    # 
    # Based on GeneType {'protein_coding', rRNA', 'tRNA', etc}, we can extract all annotated rRNA genes that are then depleted by applying `rRNAdepletion(df,rRNAs)`.

    # Load ggf3 file and filter for gene entries
    gff3 = pd.read_csv(gffPath, sep='\\t', header = None, skiprows = findRow(gffPath))
    gff3 = gff3.loc[gff3.iloc[:,2] == 'gene']

    # Format new columns GeneID, GeneType and GeneSymbol by extracting substrings of attributes column
    gff3['ID'] = pd.DataFrame(gff3.iloc[:,8].str.split('ID=', expand = True)).iloc[:,1].str.split(';', expand = True).iloc[:,0]
    gff3['GeneType'] = pd.DataFrame(gff3.iloc[:,8].str.split('gene_biotype=', expand = True)).iloc[:,1].str.split(';', expand = True).iloc[:,0]
    gff3['Symbol'] = pd.DataFrame(gff3.iloc[:,8].str.split('gene=', expand = True)).iloc[:,1].str.split(';', expand = True).iloc[:,0]

    # Add entity host and phage column based on seqid column
    # entity = {'NC_000913.3' : 'host', 'NC_000866.4' : 'phage'} # this needs to be adapted depending on phage-host interaction

    # automation assuming there are more host genes than phage genes in dual genome
    entity = findEntity(gff3)
    gff3['Entity'] = gff3[0].map(entity)

    # Index dataframe with GeneID (important for rRNAdepeltion)
    gff3.index = gff3['ID']

    # Extract list with rRNA entries based on GeneType column
    rRNAs = gff3.loc[gff3['GeneType'] == 'rRNA', 'ID'].tolist()

    # Perform rRNA depletion
    df_norRNAs = rRNAdepletion(df,rRNAs)

    # Verify rRNA depletion
    # Maybe throw error?
    # print(f"Genes before rRNA depletion: {df.shape[0]}")
    # print(f"Genes after  rRNA depletion: {df_norRNAs.shape[0]}")

    # ### 2.3 TPM normalization

    # Now we can tpm normalize the raw counts table using `TPM(df_norRNAs, df_initial, 0.5)` specifying a pseudocount of `0.5`. We further add `Entity` (Host or Phage) and `Symbol` (GeneSymbol) to the `tpms` and `df_norRNAs` dataframe.

    # TPM normalization of raw counts
    tpms = TPM(df_norRNAs, df_initial, 0.5)

    # Add Entity and Symbol columns to df_norRNAs
    tpms['Entity'] = gff3.loc[sorted(tpms.index.to_list()), 'Entity'].astype(str)
    tpms['Symbol'] = gff3.loc[sorted(tpms.index.to_list()), 'Symbol'].astype(str)
    df_norRNAs[['Entity', 'Symbol']] = tpms[['Entity', 'Symbol']]

    # Next, we use `fillSymbols(tpms)` to assign gene IDs to entries where gene symbols are missing (i.e., `Symbol = None`). Then, `make_unique_with_index(tpms)` is applied to ensure GeneSymbols are unique. If not, corresponding GeneIDs are concatenated with GeneSymbols.


    # GeneSymbols are added to tpm table
    tpms = fillSymbols(tpms)
    tpms = make_unique_with_index(tpms)

    # ### 2.4 logTPM normalization

    # TPM values are further log-transformed using `logNorm(tpms.iloc[:, :-2])` (raw TPM matrix as input). Entity and Symbols are added afterwards.

    # Log transformation of tpm values
    logTPMs = logNorm(tpms.iloc[:, :-2])
    logTPMs = logTPMs.join(tpms.iloc[:, -2:])

    # ## 3 Filter samples, if necessary
    # 
    # After tpm normalization, we now evaluate data quality by performing PCA anaylsis on our replicates (tpms and logTPMs) using the `txPCA()` function, which again takes the raw Data matrix without `Entity` and `Symbol` as input. Potential outliers can be identified here and filtered out for further data analysis.

    # This needs to be adapted depending on the dataset

    # columnOrder = ['0_R1', '0_R2', '0_R3', 
    #                '1_R1', '1_R2', '1_R3', 
    #                '4_R1', '4_R2', '4_R3', 
    #                '7_R1', '7_R2', '7_R3', 
    #                '20_R1', '20_R2', '20_R3']

    columnOrder = df.sort_index(axis=1).columns.tolist()

    # PCA on tpms values
    txPCA(tpms[columnOrder])
    plt.title("PCA on TPM values")

    # PCA on logTPM values
    txPCA(logTPMs[columnOrder])
    plt.title("PCA on logTPM values")

    # For further analysis, PCA is performed on host and phage genes separately. For this, we filter our dataframes for the corresponding genes.

    tpmHost = tpms[tpms['Entity'] == 'host']
    tpmPhage = tpms[tpms['Entity'] == 'phage']
    logTPMsHost = logTPMs[logTPMs['Entity'] == 'host']
    logTPMsPhage = logTPMs[logTPMs['Entity'] == 'phage']

    # tpmHost
    txPCA(tpmHost[columnOrder])
    plt.title("PCA on TPM_Host values")

    # tpmPhage
    txPCA(tpmPhage[columnOrder])
    plt.title("PCA on TPM_Phage values")

    # logTPM Host
    txPCA(logTPMsHost[columnOrder])
    plt.title("PCA on logTPM_Host values")

    # logTPM Phage
    txPCA(logTPMsPhage[columnOrder])
    plt.title("PCA on logTPM_Phage values")

    # Update all tables by excluding 4_R1
    # updatedOrder = ['0_R1', '0_R2', '0_R3', '1_R1', '1_R2', '1_R3', '4_R2', '4_R3', '7_R1', '7_R2', '7_R3', '20_R1', '20_R2', '20_R3']
    # updatedOrder = df.drop(['4_R1'], axis = 1).sort_index(axis=1).columns.tolist()

    # Or if no samples need to be excluded

    # Keep tables as is if no exclusion necessary 
    updatedOrder = columnOrder

    # Updated tpms
    txPCA(tpms[updatedOrder])
    plt.title("PCA on TPM values after excluding sample 4_R1")

    # Updated logTPMs
    txPCA(logTPMs[updatedOrder])
    plt.title("PCA on logTPM values after excluding sample 4_R1")

    # Updated tpmHost
    txPCA(tpmHost[updatedOrder])
    plt.title("PCA on TPM_Host values after excluding sample 4_R1")

    # tpmPhage
    txPCA(tpmPhage[updatedOrder])
    plt.title("PCA on TPM_Phage values after excluding sample 4_R1")

    # logTPM Host
    txPCA(logTPMsHost[updatedOrder])
    plt.title("PCA on logTPM_Host values after excluding sample 4_R1")

    # logTPM Phage
    txPCA(logTPMsPhage[updatedOrder])
    plt.title("PCA on logTPM Phage values after excluding sample 4_R1")

    # After PCA analysis on data with excluded 4_R1, we decided to remove this datapoint from further analysis.

    # ## 4. Final grouping
    # 
    # We will now summarize time points with mean and standard deviation for TPM-normalized data. First, we update the tpms and raw counts dataframe by dropping sample 4_R1, as well as adding the Entity and Symbol column to df_norRNAs.

    # Updated tpms dataframe 
    tpmsUpdated = tpms #.drop('4_R1', axis = 1)

    # Updated raw counts dataframe
    df_norRNAsUpdated = df_norRNAs #.drop('4_R1', axis = 1)

    # Next, we generate dataframes containing mean and standard deviation of gene expression values at specific timepoints by applying `getMeanSD` which again takes raw tpms matrix (without `Entity` and `Symbol`) as input. After creation, we add `Entity` and `Symbol` again.

    # Create dataframes for tpm means and standard deviations
    TPMmeans, TPMsds = getMeanSD(tpmsUpdated[updatedOrder])

    # Add Entity and Symbol to TPMmeans again
    # TPMmeans = TPMmeans[['0', '1', '4', '7', '20']] # adapt to dataset
    TPMmeans = TPMmeans.sort_index(axis=1)
    TPMmeans[['Entity', 'Symbol']] = tpms[['Entity', 'Symbol']]

    # Add Entity and Symbol to TPMsds again
    # TPMsds = TPMsds[['0', '1', '4', '7', '20']] # adapt to dataset
    TPMsds = TPMsds.sort_index(axis=1)
    TPMsds[['Entity', 'Symbol']] = tpms[['Entity', 'Symbol']]

    # We also calculate relative expression values for each gene, using TPM-normalized means to scale to highest expression per gene across time points, by running `proportionalExp` on the raw TPMmeans matrix. Afterwards, we again add `Entity` and `Symbol` to the dataframe.

    # Create dataframe with relative expression values compared to highest expression timepoint
    # TPMmeansCopy = TPMmeans[['0','1','4','7','20']].copy() # adapt to dataset
    TPMmeansCopy = TPMmeans[findTimes(TPMmeans)].copy()
    propExp = proportionalExp(TPMmeansCopy)
    propExp[['Entity', 'Symbol']] = TPMmeans[['Entity', 'Symbol']]
    propExp.head(2)

    # ## 5. Phage gene classification
    # 
    # Next, we classify phage genes into early, middle, or late categories based on their temporal expression profiles throughout the infection. For this, we use two functions — `classLabelThreshold()` and `classLabelMax()` — which assign genes to classes based on:
    # 
    # - whether their expression at a given timepoint exceeds 20% of their maximum expression (classLabelThreshold), or
    # - the timepoint at which maximum expression occurs (classLabelMax).
    # 
    # Timepoint boundaries for classification can chosen based on the given data as well as already known criteria. They are hardcoded into the classification functions:
    # 
    # - Early genes: Meet the classification criteria at timepoints 1 or 4
    # - Middle genes: Meet the classification criteria at timepoint 7
    # - Late genes: Meet the classification criteria at timepoint 20
    # 
    # One may need to add more conditions, if more timepoints are available.

    # Classify phage genes according to set timepoints in the respective functions
    TPMmeans = classLabelThreshold(TPMmeans,5,10)
    TPMmeans = classLabelMax(TPMmeans,5,10)

    TPMmeans[TPMmeans["Entity"] == "phage"].head(2)

    # We also add determined classes to other dfs.

    TPMsds[['ClassThreshold', 'ClassMax']] = TPMmeans[['ClassThreshold', 'ClassMax']]
    tpmsUpdated[['ClassThreshold', 'ClassMax']] = TPMmeans[['ClassThreshold', 'ClassMax']]
    propExp[['ClassThreshold', 'ClassMax']] = TPMmeans[['ClassThreshold', 'ClassMax']]

    # ## 6. Add variance to all dataframes

    # Lastly, we will assess the variability of each gene based on its gene expressionl. To do so, we calculate the stabilized variance of tpm expression values (variance / mean) for each gene across all timepoints and replicates using the function `stabilizedVariance`.

    # Add stabilized variance of genes over timepoints to tpms dataframe
    tpmsUpdated = stabilizedVariance(tpmsUpdated)
    tpmsUpdated.head()

    # Calculated variances are further added to other dataframes.

    # Add variance to other dataframes
    TPMmeans['Variance'] = tpmsUpdated['Variance']
    TPMsds['Variance'] = tpmsUpdated['Variance']
    propExp['Variance'] = tpmsUpdated['Variance']

    # ## 7. Write data to output

    # Full TPM table
    tpmsUpdated.to_csv('Wolfram-Schauerte_full_TPM.tsv', sep = '\\t')
    # Full raw_counts table
    df_norRNAsUpdated.to_csv('Wolfram-Schauerte_full_raw_counts.tsv', sep = '\\t')
    # Summarized (time point means) TPM table
    TPMmeans.to_csv('Wolfram-Schauerte_TPM_means.tsv', sep = '\\t')
    # Summarized (time point) TPM standard deviation
    TPMsds.to_csv('Wolfram-Schauerte_TPM_std.tsv', sep = '\\t')
    # Proportional expression per gene and time point
    propExp.to_csv('Wolfram-Schauerte_fractional_expression.tsv', sep = '\\t')
    """
}