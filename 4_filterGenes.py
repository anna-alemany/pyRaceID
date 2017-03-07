#### Dependencies ####
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
import argparse as args
import pickle 

#### Input parameters/options ####
parser = args.ArgumentParser(description = 'Filtering raw couttables using QC, run after 2_filterCells')
# files
parser.add_argument('--picklein', help = 'input pickle file')
parser.add_argument('--out', help = 'output file (without extension)', default = 'None')

# filter numbers
parser.add_argument('--mincells', type = float, help = 'min. number of cells where gene is present (def=2)', default = 2)
parser.add_argument('--spikein', help = 'when called, spikein genes are removed from dataset', action = 'store_true')
parser.add_argument('--chrM', help = 'when called, mitochondrial genes are removed', action = 'store_true')
parser.add_argument('--minexpr', type = float, help = 'min. total expression per gene (def = 5)', default = 5)
parser.add_argument('--maxexpr', type = float, help = 'max. total expression per gene (def = infty)', default = np.infty)

# conversion to script parameters
args = parser.parse_args()

inputpickle = args.picklein
if not os.path.isfile(inputpickle):
    sys.exit('input pickle not found')

outputfile = args.out

if outputfile == 'None':
    outputfile = inputpickle[:-7] + '_filteredGenes'
    foutputfile = inputpickle[:-7]

mincells = args.mincells
spikein = args.spikein
chrm = args.chrM
minexpr = args.minexpr
maxexpr = args.maxexpr

#### Filtering ####
pdata = pickle.load(open(inputpickle))

ndata = pdata['normDat']

fdata = ndata[(ndata >= minexpr).sum(axis = 1) >= mincells]
fdata = fdata[fdata.max(axis=1) < maxexpr]

if spikein:
    fdata = fdata.loc[['ERCC' not in idx for idx in fdata]]

if chrm:
    fdata = fdata.loc[['chrM' not in idx for idx in fdata]]
    fdata = fdata.loc[['mt-' not in idx for idx in fdata]]

fdata.to_csv(outputfile + '.txt', sep = '\t')

pdata['filteredDat'] = fdata
pdata['filterGenesParametres'] = {'mincells': mincells,
        'spikein': spikein,
        'chrM': chrm,
        'minexpr': minexpr, 
        'maxexpr': maxexpr}

pickle.dump(pdata, open(foutputfile + '.pickle', 'w'))




