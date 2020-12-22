in_dir='./data/transformed/1_reformatted'
out_dir='./data/transformed/2_clean'
misc_dir='./data/misc/'

mkdir -p $out_dir
mkdir -p $misc_dir

#--------------------------------------------------#
# Creating ontology mappings

## sample.combined_exc_4_species_integration.h5ad
file='sample.combined_exc_4_species_integration.h5ad'
python3 ./scripts/create_ontology_lookup.py $in_dir/$file ./data/misc/$file CL:0000679

# sample.combined_exc_integration.h5ad
file='sample.combined_exc_integration.h5ad'
python3 ./scripts/create_ontology_lookup.py $in_dir/$file ./data/misc/$file CL:0000679

# sample.combined_glia_integration.h5ad
file='sample.combined_glia_integration.h5ad'
python3 ./scripts/create_ontology_lookup.py $in_dir/$file ./data/misc/$file non-neuronal

# sample.combined_glia_integration.h5ad
file='sample.combined_inh_integration.h5ad'
python3 ./scripts/create_ontology_lookup.py $in_dir/$file ./data/misc/$file CL:0000617

#--------------------------------------------------#
# Appending BICCN ontology and rename columns
for i in $in_dir/*h5ad
do
    base=$(basename $i)
    out_file=$out_dir/$base
    if [ ! -f $out_file ]
    then
        python3 scripts/append_biccn_ontology.py $i ./data/misc/ontology_biccn.txt $out_dir/temp_$base
        python3 scripts/append_sex.py $out_dir/temp_$base ./data/misc/sex_mappings.tsv $out_dir/temp2_$base
        python3 scripts/append_organism.py $out_dir/temp2_$base ./data/misc/${base}_ontology_lookup_organism.tsv $out_dir/temp3_$base
        python3 scripts/rename_columns_anndata.py cluster_label:BICCN_cluster_label,subclass_label:BICCN_subclass_label $out_dir/temp3_$base $out_file && rm $out_dir/temp*_$base
    fi
done
