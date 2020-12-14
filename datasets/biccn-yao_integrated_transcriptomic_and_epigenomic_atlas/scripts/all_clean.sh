in_dir='./data/transformed/1_reformatted'
out_dir='./data/transformed/2_cleaned'

mkdir -p $out_dir

for i in $in_dir/*_filtered*
do
    
    base=$(basename $i)
    out_file=$out_dir/$base
    
    if [ ! -f $out_file ]
    then
        echo "Submmitted $base"
        python3 scripts/filter_nan.py $i class_label $out_file &
    fi
done
