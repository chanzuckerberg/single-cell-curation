out_dir='./data/transformed/4_remixed'
in_dir='./data/transformed/3_prepared'
ontology_table='./data/misc/ontology_lookup_cell_type.tsv'
yaml_out='./schema_config_files'

mkdir -p $out_dir
mkdir -p $yaml_out

# Create config files
Rscript scripts/create_yaml.R $ontology_table "10X_cells_v2_AIBS" "EFO:0009899" "Male" "Female" $yaml_out/10X_cells_v2_AIBS.yml
Rscript scripts/create_yaml.R $ontology_table "10X_cells_v3_AIBS" "EFO:0009922" "Male" "Female" $yaml_out/10X_cells_v3_AIBS.yml
Rscript scripts/create_yaml.R $ontology_table "10X_nuclei_v3_AIBS" "EFO:0009922" "Male" "Female" $yaml_out/10X_nuclei_v3_AIBS.yml
Rscript scripts/create_yaml.R $ontology_table "10X_nuclei_v3_Broad" "EFO:0009922" "MALE" "FEMALE" $yaml_out/10X_nuclei_v3_Broad.yml
Rscript scripts/create_yaml.R $ontology_table "SMARTer_cells_MOp" "EFO:0008930" "M" "F" $yaml_out/SMARTer_cells_MOp.yml
Rscript scripts/create_yaml.R $ontology_table "SMARTer_nuclei_MOp" "EFO:0008930" "M" "F" $yaml_out/SMARTer_nuclei_MOp.yml

#
out_file=$out_dir/SMARTer_cells_MOp_filtered.h5ad
in_file=$in_dir/SMARTer_cells_MOp_filtered.h5ad
if [ ! -f $out_file ]
then
    echo "Submmitted $out_file"
    cellxgene schema apply --source-h5ad $in_file --remix-config $yaml_out/SMARTer_cells_MOp.yml --output-filename $out_file
fi


out_file=$out_dir/SMARTer_nuclei_MOp_filtered.h5ad
in_file=$in_dir/SMARTer_nuclei_MOp_filtered.h5ad
if [ ! -f $out_file ]
then
    echo "Submmitted $base"
    cellxgene schema apply --source-h5ad $in_file --remix-config $yaml_out/SMARTer_nuclei_MOp.yml --output-filename $out_file
fi

out_file=$out_dir/10X_cells_v2_AIBS_filtered.h5ad
in_file=$in_dir/10X_cells_v2_AIBS_filtered.h5ad
if [ ! -f $out_file ]
then
    echo "Submmitted $base"
    cellxgene schema apply --source-h5ad $in_file --remix-config $yaml_out/10X_cells_v2_AIBS.yml --output-filename $out_file
fi

#
out_file=$out_dir/10X_cells_v3_AIBS_filtered.h5ad
in_file=$in_dir/10X_cells_v3_AIBS_filtered.h5ad
if [ ! -f $out_file ]
then
    echo "Submmitted $base"
    cellxgene schema apply --source-h5ad $in_file --remix-config $yaml_out/10X_cells_v3_AIBS.yml --output-filename $out_file
fi

#
out_file=$out_dir/10X_nuclei_v3_AIBS_filtered.h5ad
in_file=$in_dir/10X_nuclei_v3_AIBS_filtered.h5ad
if [ ! -f $out_file ]
then
    echo "Submmitted $base"
    cellxgene schema apply --source-h5ad $in_file --remix-config $yaml_out/10X_nuclei_v3_AIBS.yml --output-filename $out_file
fi

#
out_file=$out_dir/10X_nuclei_v3_Broad_filtered.h5ad
in_file=$in_dir/10X_nuclei_v3_Broad_filtered.h5ad
if [ ! -f $out_file ]
then
    echo "Submmitted $base"
    cellxgene schema apply --source-h5ad $in_file --remix-config $yaml_out/10X_nuclei_v3_Broad.yml --output-filename $out_file
fi

