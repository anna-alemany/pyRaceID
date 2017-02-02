# pyRaceID
RaceID in python.

## Possible pipeline 

### 1. Quality Controls

Script: QC.py

Input: coutt, coutc and coutb tables. The scripts asks you for the root of the name of the tables, supposed to be the same for the three tables.

Checks total sum per cell per table (sum-b, sum-c, sum-t), total number of genes per cell (gene-t; includes MT and Spike), fraction of mitochondria genes (mtfrac-t), number of spike-in per cell (spikein-t), overexpression per cell (overexpr), fraction of SvdB stressed genes (SvdBfracsg).

Output: root-QC.txt (table: rows: parameters desribed above; columns: cells)

Data visualization: I have preference for gnuplot. <br/>
To produce violing plots of the data for a given column go into gnuplot terminal and type (more info here: http://gnuplot.sourceforge.net/demo_cvs/violinplot.html):

col = 2  <br/>
set table $violin  <br/>
pl 'root-QC.txt' us col:(1) smooth kdensity  <br/>
unset table  <br/>
pl $violin us 2:1 w l lw 2 lc 1  noti<br/>
repl $violin us (-$2):1 w l lw 2 lc 1 noti <br/>
repl 'root-QC.txt' us (0):col:(.5) with boxplot ti col  <br/>

### 2. Filtering

Based on QC one can decide which cells to filter.
There is no script, since it can be vary from set to set.

#### 3. Normalization

Methods:
- Library size:
- TMM: (edgeR package)
- DESeq: (DESeq package)
- Downsampling:
- Spikein:
- Deconvolution:



