import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


# Adding annotation to metadata, matching samples in raw counts table

def annotateData(df, sampleDict):
    dfOut = df.copy()
 
    dfOut['SampleID'] = dfOut['Run'].astype(str) + '_sorted.bam'

    dfOut.index = dfOut['Sample Name']

    dfOut['SampleNames'] = pd.Series(sampleDict)
    dfOut.index = dfOut['SampleID']

    return dfOut

# Changing column names based on annotated metadata

def changeColnames(df, meta):
    newCols = meta.loc[df.columns]['SampleNames']
    df.columns = newCols
    df = df.reindex(sorted(df.columns), axis =1)

    return df

# Perform in silico rRNA depletion

def rRNAdepletion(df, rRNAs):
    genes = list(set(df.index) - set(rRNAs))
    df = df.loc[genes,:]
    return df

# Function for conversion to TPM and addition of pseudocount

def TPM(df, meta, pse):
    """
    df: rRNA depleted input table
    meta: original dataframe including 'Length' column
    pse: pseudocount
    """
    
    lengths = meta.loc[df.index,'Length']
    tpmData = df.copy()
    tpmData = tpmData.astype('float64')  

    for i in range(0,tpmData.shape[1]):
        rpk = (tpmData.iloc[:,i]+pse)/lengths
        scalingfactor = np.sum(rpk)/1000000
        tpm = rpk/scalingfactor

        tpmData.iloc[:,i] = tpm
    
    return tpmData

# Function for conversion to log2(x+1)

def logNorm(df):
    """
    df: rRNA depleted input table
    """
    
    logData = df.copy()
    logData = np.log2(np.array(logData)+1) # add 1 to have at least always positive values
    logData = pd.DataFrame(logData, index = df.index, columns = df.columns)
    
    return logData

# Function for PCA and visualization

def txPCA(df):
    Df = df.transpose().to_numpy()
    pca = PCA(n_components=2)
    components = pca.fit_transform(Df)
    X = pd.DataFrame(pd.concat([pd.DataFrame(components, index=df.columns), pd.DataFrame(df.columns.tolist(), index=df.columns)], axis = 1, ignore_index=False))
    X.columns = ["Dim1", "Dim2", "Sample"]
    X['TimePoint'] = X['Sample'].str.split('_', expand=True).iloc[:,0]
    X['Replicate'] = X['Sample'].str.split('_', expand=True).iloc[:,1]
    ax = sns.scatterplot(x = "Dim1", y = "Dim2", data = X, hue = "TimePoint", style = "Replicate")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_xlabel("Dim 1 (explained variance " + str(round(pca.explained_variance_ratio_[0] *100, ndigits=2)) + " %)")
    ax.set_ylabel("Dim 2 (explained variance " + str(round(pca.explained_variance_ratio_[1] *100, ndigits=2)) + " %)")

# Getting mean and sd for samples of same time points

def getMeanSD(df):
    # Get samples
    sampleBase = set([x[:-1] for x in df.columns])
    samples = [x[:-1] for x in df.columns]
    
    # Define a new dataframe to store mean and sd
    newCols = [x[:-2] for x in sampleBase]
    means = df.iloc[:,:len(newCols)]
    means.columns = newCols

    newCols = [x[:-2] for x in sampleBase]
    sds = df.iloc[:,:len(newCols)]
    sds.columns = newCols

    # Loop over sample base
    for base in sampleBase:
        
        indices = [i for i in range(0,len(samples)) if base == samples[i]]
        mean = base[:-2]
        sd = base[:-2]

        means[mean] = np.mean(df.iloc[:,indices], axis = 1)
        sds[sd] = np.std(df.iloc[:,indices], axis = 1)

    return pd.DataFrame(means), pd.DataFrame(sds)

# Calculation of stabilized variance (expression value adjusted)

def stabilizedVariance(df):
    labels = list()
    
    i = 0
    while i < df.shape[0]:

        # Get array of expression values at time points
        expressions = list(df.iloc[i,0:(df.shape[1]-4)])

        # Get mean expression for the gene
        exprMean = np.mean(np.array(expressions))

        # Get the variance for the gene
        varGene = np.var(np.array(expressions))

        # Stabilized variance
        stableVarGene = varGene/exprMean

        labels.append(stableVarGene)

        i += 1

    tpmOut = df.copy()
    tpmOut['Variance'] = labels

    return tpmOut

# Use TPM-normalized means to scale to highest expression per gene across time points

def proportionalExp(df):
    normExp = df.copy()
    for i in range(normExp.shape[0]):
        maxValue = np.max(np.array(normExp.iloc[i,:]))
        normExp.iloc[i,:] = normExp.iloc[i,:]/maxValue

    return normExp

# Function to fill in missing symbols by geneid.

def fillSymbols(df):
    df_new = df.copy()
    index = df.index.to_list()
    for i in range(0,df.shape[0]):
        if (df.iloc[i,-1:].values == None):
            df_new.iloc[i,-1:] = index[i]
    return df_new

# Make gene symbols unique using index (gene IDs)

def make_unique_with_index(df):
    count_dict = {}
    unique_list = []
    
    lst = df['Symbol'].tolist()
    index = df.index.tolist()

    for i in range(0, len(lst)):
        item = lst[i]
        if item in count_dict:
            count_dict[item] += 1
            new_item = f"{index[i]}_gene_{item}"
        else:
            count_dict[item] = 0
            new_item = item
        
        unique_list.append(new_item)
    
    df_new = df
    df_new['Symbol'] = unique_list
    
    return df_new

# Make gene symbols unique

def make_unique(lst):
    count_dict = {}
    unique_list = []
    
    for item in lst:
        if item in count_dict:
            count_dict[item] += 1
            new_item = f"{item}_{count_dict[item]}"
        else:
            count_dict[item] = 0
            new_item = item
        
        unique_list.append(new_item)
    
    return unique_list

# Add a classification label based on exceeding 20 % of maximal expression
# Takes time after which genes are "middle" or "late" genes

def classLabelThreshold(tpm, middle, late):
    labels = []
    #extract all measured timings from column names
    timings = tpm.iloc[:,:-2].columns.astype(float).tolist()

    i = 0
    while i < tpm.shape[0]:
        # Get array of expression values at time points
        expressions = list(tpm.iloc[i,0:(tpm.shape[1]-2)])
        # Get maximal value for each gene across time points
        maxTPM = max(expressions)
        # Get the threshold value
        thresHold = maxTPM*0.2
        # get index of time point
        indices = [x for x in range(len(expressions)) if expressions[x] >= thresHold]
        # convert index to time point
        timePoint = timings[min(indices)]

        if timePoint == 0:
            labels.append('None')
        elif timePoint < middle:
            labels.append('early')
        elif timePoint < late:
            labels.append('middle')
        elif timePoint >= late:
            labels.append('late')
        else:
            labels.append('NaN')

        i += 1

    tpmOut = tpm.copy()
    tpmOut['ClassThreshold'] = labels
    tpmOut.loc[tpmOut['Entity'] == 'host', 'ClassThreshold'] = 'None'

    return tpmOut


# Add a classification label based on maximal expression 
def classLabelMax(tpm, middle, late):
    labels = []
    #extract all measured timings from column names
    timings = tpm.iloc[:,:-3].columns.astype(float).tolist()

    i = 0
    while i < tpm.shape[0]:
        # Get array of expression values at time points
        expressions = list(tpm.iloc[i,0:(tpm.shape[1]-3)])
        # Get maximal value for each gene across time points
        maxTPM = max(expressions)
        # Get the threshold value
        thresHold = maxTPM*0.2
        # get index of time point
        indices = [x for x in range(len(expressions)) if expressions[x] >= thresHold]
        # convert index to time point
        timePoint = timings[min(indices)]

        if timePoint == 0:
            labels.append('None')
        elif timePoint < middle:
            labels.append('early')
        elif timePoint < late:
            labels.append('middle')
        elif timePoint >= late:
            labels.append('late')
        else:
            labels.append('NaN')

        i += 1

    tpmOut = tpm.copy()
    tpmOut['ClassMax'] = labels
    tpmOut.loc[tpmOut['Entity'] == 'host', 'ClassMax'] = 'None'

    return tpmOut