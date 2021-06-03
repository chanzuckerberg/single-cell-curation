mkdir -p ./data/original/

cd ./data/original

wget https://github.com/berenslab/mini-atlas/raw/master/data/m1_patchseq_ephys_features.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/m1_patchseq_exon_counts.csv.gz
wget https://github.com/berenslab/mini-atlas/raw/master/data/m1_patchseq_intron_counts.csv.gz
wget https://github.com/berenslab/mini-atlas/raw/master/data/m1_patchseq_meta_data.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/m1_patchseq_morph_features.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/m1_patchseq_morph_zprofiles.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/gene_lengths.txt
wget https://github.com/berenslab/mini-atlas/raw/master/data/processed/export/coverage-tsnes.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/processed/export/morpho-electric-tsnes.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/processed/export/morpho-electric-tsnes.csv
wget https://github.com/berenslab/mini-atlas/raw/master/data/processed/rnaseq/10x-tsne-pvsst.pickle
wget https://github.com/berenslab/mini-atlas/raw/master/data/processed/rnaseq/10x-tsne-exc.pickle
wget https://github.com/berenslab/mini-atlas/raw/master/data/processed/rnaseq/10x-tsne-viplamp.pickle
wget -O 10X_cells_v2_AIBS.pickle.zip https://www.dropbox.com/s/tkbrn6egwq111iu/10X_cells_v2_AIBS.pickle.zip?dl=0
curl -o 10X_cells_v2_AIBS.h5ad "https://corpora-data-prod.s3.amazonaws.com/c88e2a9c-72b8-4a88-a2f6-e428eada0c86/local.h5ad?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=ASIATLYQ5N5X4DOOYDG4%2F20210427%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20210427T172520Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Security-Token=IQoJb3JpZ2luX2VjENH%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQC9hAV2Gy7OArqOw7v6HYvQ0OGSYv82NBfSsMQ6aY6xdgIhALDmsaCcu8yRWQgJ%2FKvwon%2F%2FG1%2BNskyBZ3WRe%2Fflk%2Bs7KqUDCEkQARoMMjMxNDI2ODQ2NTc1IgwFVQfAcre7eatF2gsqggMHBjZlfY2omQiKEXk84tUxv0t4F5OgyxeeMmZuhP%2BAfvif0H6TLASLUH1nwnaBxiJz8fUo3QXTrreBDOzmTumK4lofNSFSH0UZTuMEG6Ru%2BHKBtsGi93pj0hebGN7PPF6L9LpJ1Wq%2F2rUd%2FCMoBvmqEb9hAnlU%2F6MVpQ5nu%2Fr8fJ%2ByyqT8o8JLNYu0xXUrnBzZD0KYxb1PE2teLKlKZWDykHi0LfHM3V1GdGt8nZajp3MB53hdFgtAW6NDoxzE3A7a1EsMWXoyK%2Bl0EHyNAiEJ1mpMaU9JOF0lx0qb7Yl9JkyUA9H6xXgIWdvbgah9RhPZ2qzejutwZsUqHsZTu5f8yFjBBOD4YpLzDYZ7YnWjsTqglt2okJaKgQ%2Bk7vl4tfxtp%2FkbCbgJO1BsqLJIF%2FS3X4UizR%2BFY6m8nHfdiX7C3CFeAhq5ECVy8lJHjhnkH2O7rVk1mf7%2BDCMDAH55iNlNICwkK5PIdyj90%2BkhMqIFaHMOKdRKOKaee7oX6f%2BY8g2XPjCT86CEBjrqAcN7OsgZCBDmQa6t5YEgKiCsYEWyIZOMPkLi58t4NF5Q66wyv7nseYEhX7J54WbuQhvSKwAvuofBbphBamrLcQKHogfByOPhTzyLxiz%2B5KjTZXdIe%2BERiOkdn9zqNKBjhe5Xbc12sGHjK7iJZ7XWryXkpSdSlJwgJhjAWtkdCyDamWsFVhQOHWtE%2FxFdf5BQ9r5%2Bj3yB47WjnYoRoI5z9XjyocGZHuvpjwTqzw7V9K1VA8iXop8yKz1Ivag5LajLPMjvcw%2FnuM%2FI3OaO%2F9Zbmepc8GJXYqE0J4cDYzvS%2Bdb6kgktt8xLiYRsbw%3D%3D&X-Amz-Signature=562a20e7df303cd08492ae3bfe59167068111d5651e9c83ce7182d3c1cbee483"
unzip 10X_cells_v2_AIBS.pickle.zip 
rm 10X_cells_v2_AIBS.pickle.zip 
