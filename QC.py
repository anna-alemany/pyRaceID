import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv

try:
    rinput = sys.argv[1]
except:
    sys.exit('Please, give root to cout tables')

#### SvdB stressed genes ####
SvdBsg = ["Actg1__chr11","Ankrd1__chr19","Arid5a__chr1","Atf3__chr1","Atf4__chr15","Bag3__chr7","Bhlhe40__chr6","Brd2__chr17",
        "Btg1__chr10","Btg2__chr1","Ccnl1__chr3","Ccrn4l__chr3","Cebpb__chr2","Cebpd__chr16","Cebpg__chr7","Csrnp1__chr9",
        "Cxcl1__chr5","Cyr61__chr3","Dcn__chr10","Ddx3x__chrX","Ddx5__chr11","Des__chr1","Dnaja1__chr4","Dnajb1__chr8",
        "Dnajb4__chr3","Dusp1__chr17","Dusp8__chr7","Egr1__chr18","Egr2__chr10","Eif1__chr11","Eif5__chr12","Erf__chr7",
        "Errfi1__chr4","Fam132b__chr1","Fos__chr12","Fosb__chr7","Fosl2__chr5","Gadd45a__chr6","Gcc1__chr6","Gem__chr4",
        "H3f3b__chr11","Hipk3__chr2","Hsp90aa1__chr12","Hsp90ab1__chr17","Hspa1a__chr17","Hspa1b__chr17","Hspa5__chr2","Hspa8__chr9",
        "Hspb1__chr5","Hsph1__chr5","Id3__chr4","Idi1__chr13","Ier2__chr8","Ier3__chr17","Ifrd1__chr12","Il6__chr5",
        "Irf1__chr11","Irf8__chr8","Itpkc__chr7","Jun__chr4","Junb__chr8","Jund__chr8","Klf2__chr8","Klf4__chr4",
        "Klf6__chr13","Klf9__chr19","Litaf__chr16","Lmna__chr3","Maff__chr15","Mafk__chr5","Mcl1__chr3","Midn__chr10",
        "Mir22hg__chr11","Mt1__chr8","Mt2__chr8","Myadm__chr7","Myc__chr15","Myd88__chr9","Nckap5l__chr15","Ncoa7__chr10",
        "Nfkbia__chr12","Nfkbiz__chr16","Nop58__chr1","Nppc__chr1","Nr4a1__chr15","Odc1__chr12","Osgin1__chr8","Oxnad1__chr14",
        "Pcf11__chr7","Pde4b__chr4","Per1__chr11","Phlda1__chr10","Pnp__chr14","Pnrc1__chr4","Ppp1cc__chr5","Ppp1r15a__chr7",
        "Pxdc1__chr13","Rap1b__chr10","Rassf1__chr9","Rhob__chr12","Rhoh__chr5","Ripk1__chr13","Sat1__chrX","Sbno2__chr10",
        "Sdc4__chr2","Serpine1__chr5","Skil__chr3","Slc10a6__chr5","Slc38a2__chr15","Slc41a1__chr1","Socs3__chr11","Sqstm1__chr11",
        "Srf__chr17","Srsf5__chr12","Srsf7__chr17","Stat3__chr11","Tagln2__chr1","Tiparp__chr3","Tnfaip3__chr10","Tnfaip6__chr2",
        "Tpm3__chr3","Tppp3__chr8","Tra2a__chr6","Tra2b__chr16","Trib1__chr15","Tubb4b__chr2","Tubb6__chr18","Ubc__chr5",
        "Usp2__chr9","Wac__chr18","Zc3h12a__chr4","Zfand5__chr19","Zfp36__chr7","Zfp36l1__chr12","Zfp36l2__chr17","Zyx__chr6",
        "Gadd45g__chr13","Hspe1__chr1","Ier5__chr1","Kcne4__chr1"]

####

data = {'t': read_csv(rinput + '.coutt.csv', sep = '\t', index_col = False), 
        'c': read_csv(rinput + '.coutc.csv', sep = '\t', index_col = False), 
        'b': read_csv(rinput + '.coutb.csv', sep = '\t', index_col = False)}

for d in data:
    data[d] = data[d].set_index(data[d].columns[0])
    data[d].index = [i.rstrip() for i in data[d].index]

totalsum = {d: data[d].sum() for d in data}
numgenes = {d: (data[d]>0).sum() for d in data}
mtexpr = {d: data[d].loc[['__chrM' in idx for idx in data[d].index]].sum() for d in data}
mtfrac = {d: 100*mtexpr[d]/totalsum[d] for d in data}
spikein = {d: data[d].loc[['ERCC' in idx for idx in data[d].index]].sum() for d in data}
overexpr = (data['c']/data['b']).mean()
SvdBfrac = 100*data['t'].loc[[i in SvdBsg for i in data['t'].index]].sum()/totalsum['t']

df = pd.DataFrame({
    'sum-b': totalsum['b'],
    'sum-c': totalsum['c'],
    'sum-t': totalsum['t'],
    'gene-t': numgenes['t'],
    'mtfrac-t': mtfrac['t'], 
    'spikein-t': spikein['t'],
    'overexpr': overexpr,
    'SvdBfracsg': SvdBfrac
    })
df['cellID'] = df.index

print df.head()

df.to_csv(rinput + '-QC.txt', sep = '\t', index = False)

