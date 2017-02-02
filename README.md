# pyRaceID
RaceID in python.

## Possible pipeline 

### 1. Quality Controls

Input: coutt, coutc and coutb tables. The scripts asks you for the root of the name of the tables, supposed to be the same for the three tables.

Checks total sum per cell per table (sum-b, sum-c, sum-t), total number of genes per cell (gene-t; includes MT and Spike), fraction of mitochondria genes (mtfrac-t), number of spike-in per cell (spikein-t), overexpression per cell (overexpr), fraction of SvdB stressed genes (SvdBfracsg).

Output: root-QC.txt (table: rows: parameters desribed above; columns: cells)

Data visualization: I have preference for gnuplot. <br/>
To produce violing plots of the data for a given column go into gnuplot terminal and type:

col = 2
set table $violin
pl 'root-QC.txt' us col:(1) smooth kdensity
unset table
pl $violin us col:1 w l lw 2 lc 1
repl $violin us (-$col):1 w l lw 2 lc 1
repl 'root-QC.txt' us (0):2:(.5) with boxplot





