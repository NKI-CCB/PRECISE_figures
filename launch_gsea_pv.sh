#!/bin/sh

n_pv=$1
pathways="h.all c2.cp c2.cp.kegg c2.cp.reactome c2.cp.biocarta c2.cgp"

for p in $pathways
do
	cputool -c 100 -- screen -d -m sh gsea_pv_source.sh "rnaseq" $p "cell_line" $n_pv 1000
	cputool -c 100 -- screen -d -m sh gsea_pv_target.sh "rnaseq" $p "cell_line" $n_pv 1000
	cputool -c 100 -- screen -d -m sh gsea_pv_source.sh "fpkm" $p "pdx" $n_pv 1000
	cputool -c 100 -- screen -d -m sh gsea_pv_target.sh "fpkm" $p "pdx" $n_pv 1000
done 