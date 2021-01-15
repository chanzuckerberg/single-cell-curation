outdir=./data/original
mkdir -p $outdir

wget -P $outdir https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/cell_metadata.csv
wget -P $outdir https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/genes_rpkm.csv
wget -P $outdir https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/genes_counts.csv

perl -pe 's/"//g' $outdir/genes_rpkm.csv > temp && mv temp $outdir/genes_rpkm.csv
perl -pe 's/"//g' $outdir/genes_counts.csv > temp && mv temp $outdir/genes_counts.csv
