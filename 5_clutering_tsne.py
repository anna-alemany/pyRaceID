import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
import argparse as args
from sklearn.manifold import TSNE
import sklearn.cluster as clust
import scipy.spatial.distance as dist
from scipy.stats import linregress
import pickle

#### Input options/files ####
parser = args.ArgumentParser(description = 'Cluster and tsne of filtered data.')

parser.add_argument('--picklein', help = 'pickle input file, with normalzied data')
parser.add_argument('--out', help = 'output filename', default = 'None')

parser.add_argument('--metric', help = 'metric to compute distance matrix', choices = ['pearson', 'spearman', 'euclidean', 'kendall'], default = 'pearson')
parser.add_argument('--seedtsne', help = 'seed for tsne', type = int, default = 0)
parser.add_argument('--seedclust', help = 'seed for cluster', type = int, default = 10)
parser.add_argument('--cluster', help = 'cluster method (default = kmeans)', choices = ['kmeans', 'non'], default = 'kmeans')
parser.add_argument('--clusternum', help = 'max. number of clusters (in case of GAP/sat cluster number determination), or total number of clusters (otherwise) (def = 5)', type = int, default = 5)
parser.add_argument('--GAP', help = 'perform GAP statistics to determine cluster number (def = False)', action = 'store_true')
parser.add_argument('--sat', help = 'guess cluster number with saturation (def = False)', action = 'store_true')

args = parser.parse_args()

#### Convert input parameters to script parameters ####
args = parser.parse_args()

inputfile = args.picklein
output = args.out
metric = args.metric
seedtsne = args.seedtsne
seedclust = args.seedclust
clustermethod = args.cluster
clusternum = args.clusternum
findGAP = args.GAP
saturation = args.sat

if not os.path.isfile(inputfile):
    sys.exit('input pickle not found')

#### Read data ####
if not os.path.isfile(inputfile):
    sys.exit('pickle file not found')

if output == 'None':
    output = inputfile[:-7]

pdata = pickle.load(open(inputfile))

fdata = pdata['filteredDat'] # read_csv(inputfile, sep = '\t', index_col = 0)

#### Obtain distance matrix ####
if metric in ['pearson', 'spearman', 'kendall']:
    cdata = 1. - fdata.corr(metric)

elif metric in ['euclidean']:
    cdata = fdata.apply(lambda col1: fdata.apply(lambda col2: dist.euclidean(col1, col2)))

else:
    sys.exit('unknown metric')

#### Obtain tsne ####
tsne = TSNE(n_components = 2, random_state = seedtsne).fit_transform(cdata)
tsne = pd.DataFrame({'cell': cdata.columns, 'x': [tsne[i][0] for i in range(len(tsne))], 'y': [tsne[i][1] for i in range(len(tsne))]})

#### Cluster ####
def clusterfunction(method, data, clusternum, seed):
    if method == 'kmeans':
        clusters = clust.KMeans(n_clusters = clusternum, random_state = seed).fit(data)
    else:
        sys.exit('cluster method still unknown')
    return clusters

def withinClustDisp(cdata, clusters, clustnum):
    Wk = 0.
    for j in range(clustnum):
        cellList = list(cdata.columns[clusters.labels_ == j])
        Dr = 0.
        for i0 in range(len(cellList)):
            for i1 in range(i0, len(cellList)):
                Dr += cdata.loc[cellList[i0], cellList[i1]]
        Dr /= len(cellList)
        Wk += Dr
    return Wk

if not findGAP and not saturation:
    clusters = clusterfunction(clustermethod, cdata, clusternum, seedclust)
elif findGAP and not saturation:
    pass # code GAP statistics
elif saturation:
    clusters = [0 for i in range(clusternum)]
    satdf = pd.DataFrame(columns = ['k', 'sat'])
    for i in range(clusternum): 
        k = i+1;
        clusters[i] = clusterfunction(clustermethod, cdata, k, seedclust)
        Wk = withinClustDisp(cdata, clusters[i], k)
        satdf.loc[i] = [k, Wk]
    satdf['df'] = [satdf.loc[i,'sat'] - satdf.loc[min(i+1,satdf.shape[0]-1),'sat'] for i in range(satdf.shape[0])]
    fitdf = pd.DataFrame( [np.array(linregress(satdf.loc[i:,['k','sat']]))[[0,-1]] for i in range(satdf.shape[0])], columns = ['slope', 'err'])
    satdf = satdf.merge(fitdf, how = 'outer', left_index = True, right_index = True)
    satdf['slope'] = -satdf['slope']
    clusternum = int(satdf[(satdf['df'] <= satdf['slope'] + satdf['err'])]['k'].iloc[0])
    clusters = clusters[clusternum-1]
    pdata['saturation'] = satdf
    satdf.to_csv(output + '_saturation.txt', sep = '\t')

#### Heatmap ####
cellCluster = pd.DataFrame({'cell': cdata.columns, 'cluster': clusters.labels_})
sortedCells = cellCluster.sort_values(by = 'cluster')['cell']
heatmap = cdata.loc[sortedCells, sortedCells]
heatmap.columns = [''.join(['#',c]) for c in heatmap.columns]
heatmap.to_csv(output + '_heatmap.txt', sep = '\t', index = None)
heatmap.columns = [c[1:] for c in heatmap.columns]

#### tsne with cluster info ####
tsne = tsne.merge(cellCluster, how = 'outer', on = 'cell')

#tdata = fdata.transpose()
#tdata['cell'] = tdata.index
#tsne = tsne.merge(tdata, on = 'cell')

tsne.to_csv(output + '_tsne.txt', sep = '\t', index = None)

pdata['PREheatmap'] = heatmap
pdata['PREtsne'] = tsne
pdata['cluster1Parametres'] = {'metric': metric, 
        'seedtsne': seedtsne,
        'seedclust': seedclust, 
        'clustermethod': clustermethod,
        'clusternum': clusternum,
        'findGAP': findGAP, 
        'saturation': saturation, 
        'clusters': clusters}

pickle.dump(pdata, open(output + '.pickle', 'w'))

