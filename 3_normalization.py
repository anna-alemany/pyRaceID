import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
import argparse as args
import scipy.optimize
import pickle

#### Input ####
parser = args.ArgumentParser(description = 'Normalization of filtered data')
# files
parser.add_argument('--picklein', help = 'pickle file with filtered data')
parser.add_argument('--out', help = 'Root for outputfile', default = 'None')

# normalization strategy
parser.add_argument('--norm', help = 'Normalization method', choices = ['libSize', 'DSeq', 'TMM', 'SpikeIn', 'Downsample', 'Deconv'], default = 'libSize' )
parser.add_argument('--chrM', help = 'when called, mitochondrial genes are removed', action = 'store_true')
parser.add_argument('--spikein', help = 'when called, spikein genes are removed from dataset', action = 'store_true')

args = parser.parse_args()

# conversion to parameters of the script
inputdata = args.picklein
if not os.path.isfile(inputdata):
    sys.exit('pickle not found')

output = args.out
poutput = output
if output == 'None':
    output = inputdata[:-7] + '_norm'
    poutput = inputdata[:-7]

method = args.norm
spikein = args.spikein
chrm = args.chrM

#### Lets do it!
pdata = pickle.load(open(inputdata))

fdata = pdata['filterCellsDat']

if spikein:
    fdata = fdata.loc[['ERCC' not in idx for idx in fdata.index]]
if chrm:
    fdata = fdata.loc[[('chrM' not in idx or 'mt-' not in idx) for idx in fdata.index]]

if method == 'libSize': # library size normalization (standard, DEFAULT)
    ndata = (fdata/fdata.sum())*1e6 # transcripts per milion

elif method == 'Deconv': # Lun, Aaron TL, Karsten Bach, and John C. Marioni. "Pooling across cells to normalize single-cell RNA sequencing data with many zero counts." Genome biology 17.1 (2016): 75.
    numtotalcell = fdata.shape[1]
    tdata = fdata/fdata.sum()
    Ui = tdata.sum(axis = 1)
    cellList = list(fdata.sum().sort_values()[range(0,numtotalcell,2) + range(1, numtotalcell,2)[::-1]].index)*2
    w = [40, 80, 100]
    systemEqsA = pd.DataFrame(columns = cellList[:numtotalcell])
    systemEqsB = pd.DataFrame(columns = ['Vik'])
    for sizeset in w:
        for oneset in range(numtotalcell):
            cells = cellList[oneset: oneset + sizeset]
            Vik = tdata[cells].sum(axis = 1)
            Rik = (Vik/Ui).median()
            i = len(systemEqsB.index)
            systemEqsA.loc[i,cells] = 1
            systemEqsB.loc[i,'Vik'] = Rik
    systemEqsA = systemEqsA.fillna(0)
    x = scipy.optimize.nnls(np.array(systemEqsA.fillna(0)), np.array(systemEqsB['Vik']))
    s = pd.DataFrame({'Ot-1': x[0]}, index = cellList[:numtotalcell])
    s['t'] = fdata.sum()
    s['O'] = s['t']*s['Ot-1']
    ndata = fdata/s['O'] 
        

elif method == 'DS':
    sfj = (fdata.transpose()/np.exp(np.log(fdata[(rdata>0)]).mean(axis=1))).transpose()
    sfj = sfj[sfj>0].median()
    qi = (fdata/sfj).mean(axis = 1)
    wi = (fdata/sfj).var(axis = 1)
    zi = qi*(1/sfj).mean()

    pass

elif method == 'TMM':
    pass

elif method == 'Downsample':
    pass

elif method == 'Spikein':
    pass

ndata.to_csv(output + '.txt', sep = '\t')

pdata['normDat'] = ndata
pdata['normParameters'] = {
        'method': method, 
        'chrM': chrm, 
        'spikein': spikein}

pickle.dump(pdata, open(poutput + '.pickle', 'w'))

