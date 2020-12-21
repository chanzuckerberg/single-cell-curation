in_dir='./data/transformed/2_cleaned'
out_dir='./data/transformed/2_1_cleaned_ontolgy'
ontology_table='./data/misc/ontology_biccn.tsv'
ontology_table_ebi='./data/misc/ontology_lookup_cell_type.tsv'

mkdir -p $out_dir
mkdir -p $(dirname ontology_table)

# First create ontology table
python3 scripts/create_biccn_ontology_table.py $in_dir $ontology_table
python3 scripts/create_ontology_lookup_cell_type.py $in_dir $ontology_table_ebi


for i in $in_dir/*
do
    
    base=$(basename $i)
    out_file=$out_dir/$base
    
    if [ ! -f $out_file ]
    then
        echo "Submmitted $base"
        python3 scripts/append_biccn_ontology.py $i $ontology_table $out_file &
    fi
done
