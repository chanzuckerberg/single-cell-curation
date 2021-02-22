outdir=./data/original
mkdir -p $outdir

wget -O $outdir/cell_metadata.csv https://www.dropbox.com/s/hv0xrmanzananob/cell_metadata.csv?dl=1
wget -P $outdir ftp://download.brainimagelibrary.org:8811/02/26/02265ddb0dae51de/dataset_metadata.xlsx
wget -P $outdir ftp://download.brainimagelibrary.org:8811/02/26/02265ddb0dae51de/processed_data/counts.h5ad 
wget -O $outdir/umap_embedding.csv https://www.dropbox.com/s/tk9e1xe2xft5s2e/umap_embedding.csv?dl=1
