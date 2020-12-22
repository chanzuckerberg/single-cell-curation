out_dir='./data/transformed/3_remixed'
in_dir='./data/transformed/2_clean'
ontology_dir='./data/misc/'
yaml_out='./schema_config_files'

mkdir -p $out_dir
mkdir -p $yaml_out

# Create config files
Rscript scripts/create_yaml.R $ontology_dir/sample.combined_exc_4_species_integration.h5ad_ontology_lookup_cell_type.tsv "4-species integration excitory neurons" $yaml_out/sample.combined_exc_4_species_integration.h5ad.yml

Rscript scripts/create_yaml.R $ontology_dir/sample.combined_exc_integration.h5ad_ontology_lookup_cell_type.tsv "3-species integration excitory neurons" $yaml_out/sample.combined_exc_integration.h5ad.yml

Rscript scripts/create_yaml.R $ontology_dir/sample.combined_glia_integration.h5ad_ontology_lookup_cell_type.tsv "3-species integration non-nuerons" $yaml_out/sample.combined_glia_integration.h5ad.yml

Rscript scripts/create_yaml.R $ontology_dir/sample.combined_inh_integration.h5ad_ontology_lookup_cell_type.tsv "3-species integration inhibitory neurons" $yaml_out/sample.combined_inh_integration.h5ad.yml

for i in $in_dir/*h5ad
do
    base=$(basename $i)
    out_file=$out_dir/$base
    if [ ! -f $out_file ]
    then
        echo "Submmitted $out_file"
        cellxgene schema apply --source-h5ad $i --remix-config $yaml_out/${base}.yml --output-filename $out_file
    fi
done
