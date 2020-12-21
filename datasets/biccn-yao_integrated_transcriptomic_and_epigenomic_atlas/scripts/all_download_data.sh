mkdir -p ./data/original/10X_nuclei_v3_AIBS
mkdir -p ./data/original/10X_cells_v3_AIBS
mkdir -p ./data/original/10X_cells_v2_AIBS
mkdir -p ./data/original/10X_nuclei_v2_AIBS
mkdir -p ./data/original/SMARTer_cells_MOp
mkdir -p ./data/original/SMARTer_nuclei_MOp

# 10X_nuclei_v3_AIBS
wget -nc -P ./data/original/10X_cells_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/matrix.mtx.gz
wget -nc -P ./data/original/10X_cells_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/sample_metadata.csv
wget -nc -P ./data/original/10X_cells_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/features.tsv.gz
wget -nc -P ./data/original/10X_cells_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/barcode.tsv

# 10X_nuclei_v3_AIBS
wget -nc -P ./data/original/10X_nuclei_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/matrix.mtx.gz
wget -nc -P ./data/original/10X_nuclei_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/sample_metadata.csv
wget -nc -P ./data/original/10X_nuclei_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/features.tsv.gz
wget -nc -P ./data/original/10X_nuclei_v3_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/barcode.tsv

# 10X_cells_v2_AIBS
wget -nc -P ./data/original/10X_cells_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/matrix.mtx.gz
wget -nc -P ./data/original/10X_cells_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/sample_metadata.csv
wget -nc -P ./data/original/10X_cells_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/features.tsv.gz
wget -nc -P ./data/original/10X_cells_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/barcode.tsv

# 10X_nuclei_v3_Broad
wget -nc -P ./data/original/10X_nuclei_v3_Broad http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/matrix.mtx.gz
wget -nc -P ./data/original/10X_nuclei_v3_Broad http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/sample_metadata.csv
wget -nc -P ./data/original/10X_nuclei_v3_Broad http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/features.csv
wget -nc -P ./data/original/10X_nuclei_v3_Broad http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/barcodes.csv

tail -n +2 ./data/original/10X_nuclei_v3_Broad/features.csv | awk -v OFS="\t" '{print "NA\t" $1 "\tgene expression"}' > ./data/original/10X_nuclei_v3_Broad/features.tsv
gzip -f ./data/original/10X_nuclei_v3_Broad/features.tsv

awk -v OFS="\t" -v x=0 '{print ("\"" x "\"" "," "\"" $1 "\""); x=x+1}' ./data/original/10X_nuclei_v3_Broad/barcodes.csv > ./data/original/10X_nuclei_v3_Broad/barcode.tsv



# 10X_nuclei_v2_AIBS
wget -nc -P ./data/original/10X_nuclei_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/matrix.mtx.gz
wget -nc -P ./data/original/10X_nuclei_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/sample_metadata.csv
wget -nc -P ./data/original/10X_nuclei_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/features.tsv.gz
wget -nc -P ./data/original/10X_nuclei_v2_AIBS http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/barcode.tsv

# SMARTer_cells_MOp
wget -nc -P ./data/original/SMARTer_cells_MOp http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/exon.counts.csv.gz
wget -nc -P ./data/original/SMARTer_cells_MOp http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/sample_metadata.csv.gz

# SMARTer_nuclei_MOp
wget -nc -P ./data/original/SMARTer_nuclei_MOp http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/exon.counts.csv.gz
wget -nc -P ./data/original/SMARTer_nuclei_MOp http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/sample_metadata.csv.gz
