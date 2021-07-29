# This script downloads 10X data and extracts only the gtf files

get_gtf () {

  local url="$1"
  local base_dir="$2"
  local output_file="$3"

  # Download data and copies data
  curl -o 10Xdata.tar.gz $url
  tar -xvzf 10Xdata.tar.gz

  # Rename base dir because it changes from human to mouse
  mv $base_dir 10Xdata
  cp ./10Xdata/genes/genes.gtf $output_file
  rm -r 10Xdata.tar.gz 10Xdata

}

get_gtf \
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz \
refdata-gex-mm10-2020-A \
mus_musculus.gtf

get_gtf \
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz \
refdata-gex-GRCh38-2020-A \
homo_sapiens.gtf
