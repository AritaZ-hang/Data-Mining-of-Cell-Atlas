# a shell script for pySCENIC
# for more information, please see the pySCENIC tutorial

source ~/anaconda3/bin/activate pyscenic

pyscenic grn --num_workers 96 -o `Yourpath`/adjacencies.tsv    `Yourpath`/ex_matrix.csv `Yourpath`/pyscenic_database/mm_tfs.txt
		
pyscenic ctx `Yourpath`/adjacencies.tsv `Yourpath`/pyscenic_database/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather --annotations_fname `Yourpath`/pyscenic_database/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname `Yourpath`/ex_matrix.csv --mode "dask_multiprocessing" --output `Yourpath`/regulons.csv --num_workers 96

pyscenic aucell `Yourpath`/ex_matrix.csv `Yourpath`/regulons.csv -o `Yourpath`/auc_mtx.csv --num_workers 96
