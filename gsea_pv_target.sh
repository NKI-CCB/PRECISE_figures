#!/bin/sh 
data_type="$1"
pathway="$2"
source_type="$3"
n_pv=$4

if [ "$source_type" = "cell_line" ]
then
	source_name='cl'
else
	source_name='pdx'
fi

if [ ! -z $5 ] 
then 
    n_perm=$5
else
    n_perm=1000
fi

for i in $(seq 68 $n_pv)
do
	echo $i
	java -cp gsea-3.0.jar \
		-Xmx2048m xtools.gsea.Gsea \
		-res 'expression_tumors_'$data_type'.txt' \
		-cls 'scores_pv_'$data_type'_target_'$n_pv'.cls#Factor_'$i \
		-gmx 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/'$pathway'.v6.2.symbols.gmt' \
		-collapse true \
		-mode Max_probe \
		-norm meandiv \
		-nperm $n_perm \
		-permute phenotype \
		-rnd_type no_balance \
		-scoring_scheme weighted \
		-rpt_label tumors_factor_$i \
		-metric Pearson \
		-sort abs \
		-order descending \
		-chip gseaftp.broadinstitute.org://pub/gsea/annotations/ENSEMBL_human_gene.chip \
		-create_gcts false \
		-create_svgs false \
		-include_only_symbols true \
		-make_sets true \
		-median false \
		-num 100 \
		-plot_top_x 20 \
		-rnd_seed timestamp \
		-save_rnd_lists false \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-out './output/'$pathway'_breast_all_'$data_type'_tumor_'$source_name'_'$n_pv \
		-gui false
done