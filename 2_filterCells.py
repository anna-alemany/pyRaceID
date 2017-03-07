"""Filtering data
Input files should be coutt.csv tables inside a 'cout' folder.
The output file will be placed on a filtered folder, together with a LOG file"""

#### Dependencies ####
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
import glob
import pickle
import argparse as args

#### Input parameters/options ####
parser = args.ArgumentParser(description = 'Filtering raw couttables using QC')
# files
parser.add_argument('--picklein', help = 'pickle file with rawdata and QC parameters, determined previously with 1_QC.py script')
parser.add_argument('--out', help = 'output file (without file extension). If not determined, picklein is taken', default = 'None')

# filter numbers
parser.add_argument('--mint', type = float, help = 'minimum numbers of transcripts per cell (def=0)', default = 0)
parser.add_argument('--maxt', type = float, help = 'maximum number of transcripts per cell (def=infty)', default = np.infty)
parser.add_argument('--minmito', type = float, help = 'minimum fraction of mitochondria per cell (def=0)', default = 0)
parser.add_argument('--maxmito', type = float, help = 'maximum fraction of mitochondria per cell (def=infty)', default = np.infty)
parser.add_argument('--minsvdb', type = float, help = 'minimum fraction of svdb stressed genes (def=0)', default = 0)
parser.add_argument('--maxsvdb', type = float, help = 'maximum fraction of svdb stressed genes (def=infty)', default = np.infty)

# conversion to script parameters
args = parser.parse_args()

inputpickle = args.picklein
outputfile = args.out

mintranscript = args.mint
maxtranscript = args.maxt
minmitofrac  = args.minmito
maxmitofrac = args.maxmito
minsvdbfrac = args.minsvdb
maxsvdbfrac = args.maxsvdb

#### Filtering ####
if not os.path.isfile(inputpickle):
    sys.exit('ERR001: input pickle not found')

if outputfile == 'None':
    outputfile = inputpickle[:-7] 

pdata = pickle.load(open(inputpickle))

rdata = pdata['rawdata']['t'] 
info = pdata['row-QC']

info['keep'] = (info['sum-t'] >= mintranscript) & (info['sum-t'] <=  maxtranscript) & (info['mtfrac-t'] >= minmitofrac) & (info['mtfrac-t'] <= maxmitofrac) & (info['SvdBfracsg'] >= minsvdbfrac) & (info['mtfrac-t'] <= maxsvdbfrac)

fdata = rdata[info.loc[info['keep']].index]

fdata.to_csv(outputfile + '.txt' , sep = '\t')

pdata['filterCellsDat'] = fdata
pdata['filterCellsParametres'] = {'mintrans': mintranscript,
        'maxtrans': maxtranscript, 
        'minmitofrac': minmitofrac,
        'maxmitofrac': maxmitofrac,
        'minsvdbfrac': minsvdbfrac,
        'maxsvdbfrac': maxsvdbfrac, 
        'filteredCells': list(info.loc[info['keep'] == False].index)}

pickle.dump(pdata, open(outputfile + '.pickle', 'w'))



