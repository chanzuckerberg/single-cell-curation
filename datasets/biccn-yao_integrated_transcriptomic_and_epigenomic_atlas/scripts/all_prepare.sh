in_dir='./data/transformed/2_1_cleaned_ontolgy'
umap_file='./data/original/zizhen_files/umap.2d.all2.csv'
out_dir='./data/transformed/3_prepared'

mkdir -p $out_dir

for i in $in_dir/*
do
    base=$(basename $i)
    out_file=$out_dir/$base
    
    if [ ! -f $out_file ]
    then
        echo "Submmitted $base"
        #cellxgene prepare --skip-qc -r zheng17 -o $out_file $i &
        #cellxgene prepare --sparse --skip-qc  -o $out_file $i &
        python3 scripts/process_normalize_umap.py $i $umap_file $out_file &
    fi
done
