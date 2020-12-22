in_dir='./data/original/'
out_dir='./data/transformed/1_reformatted'

mkdir -p $out_dir


for i in $in_dir/*RDS
do
    
    base=$(basename $i)
    base="${base%.*}"
    out_file=$out_dir/${base}.h5ad
    
    if [ ! -f $out_file ]
    then
        echo "Submmitted $base"
        Rscript scripts/reformat_seurat.R 'curation' $i $out_file &
    fi
done
