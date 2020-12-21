mkdir -p ./data/transformed/1_reformatted

#python3 scripts/reformat_add_metadata_10X.py ./data/original/10X_cells_v2_AIBS/matrix.mtx.gz ./data/original/10X_cells_v2_AIBS/sample_metadata.csv ./data/original/10X_cells_v2_AIBS/barcode.tsv ./data/original/10X_cells_v2_AIBS/features.tsv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4c_scRNA_10X_v2_A_metadata.tsv.gz ./data/transformed/1_reformatted/10X_cells_v2_AIBS &

python3 scripts/reformat_add_metadata_10X.py ./data/original/10X_cells_v3_AIBS/matrix.mtx.gz ./data/original/10X_cells_v3_AIBS/sample_metadata.csv ./data/original/10X_cells_v3_AIBS/barcode.tsv ./data/original/10X_cells_v3_AIBS/features.tsv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4b_scRNA_10X_v3_A_metadata.tsv.gz ./data/transformed/1_reformatted/10X_cells_v3_AIBS &

python3 scripts/reformat_add_metadata_10X.py ./data/original/10X_nuclei_v3_AIBS/matrix.mtx.gz ./data/original/10X_nuclei_v3_AIBS/sample_metadata.csv ./data/original/10X_nuclei_v3_AIBS/barcode.tsv ./data/original/10X_nuclei_v3_AIBS/features.tsv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4f_snRNA_10X_v3_A_metadata.tsv.gz ./data/transformed/1_reformatted/10X_nuclei_v3_AIBS &

#python3 scripts/reformat_add_metadata_10X_broad.py ./data/original/10X_nuclei_v3_Broad/matrix.mtx.gz ./data/original/10X_nuclei_v3_Broad/sample_metadata.csv ./data/original/10X_nuclei_v3_Broad/barcode.tsv ./data/original/10X_nuclei_v3_Broad/features.tsv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4e_snRNA_10X_v3_B_metadata.tsv.gz ./data/transformed/1_reformatted/10X_nuclei_v3_Broad &

#python3 scripts/reformat_add_metadata_smart.py ./data/original/SMARTer_cells_MOp/exon.counts.csv.gz ./data/original/SMARTer_cells_MOp/sample_metadata.csv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4a_scRNA_SMART_metadata.tsv.gz ./data/transformed/1_reformatted/SMARTer_cells_MOp &

#python3 scripts/reformat_add_metadata_smart.py ./data/original/SMARTer_nuclei_MOp/exon.counts.csv.gz ./data/original/SMARTer_nuclei_MOp/sample_metadata.csv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4d_snRNA_SMART_metadata.tsv.gz ./data/transformed/1_reformatted/SMARTer_nuclei_MOp &
