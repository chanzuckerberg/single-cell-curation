# This script downloads gtf files for hunan, mouse, and sars_cvo_2. Also ERCC ids

mouse_10X_url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
human_10X_url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
covid_ftp="https://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz"

mouse_base_dir="refdata-gex-mm10-2020-A"
human_base_dir="refdata-gex-GRCh38-2020-A"

# 10X
get_gtf_10X () {

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

echo Downloading mouse GTF
get_gtf_10X \
$mouse_10X_url \
$mouse_base_dir \
mus_musculus.gtf.gz

echo Downloading human GTF
get_gtf_10X \
$human_10X_url \
$human_base_dir \
homo_sapiens.gtf.gz

# COVID-19
echo Downloading sars_cov_2 GTF
curl -o sars_cov_2.gtf.gz ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz

# ERCC
echo Downloading ERCC gene ids
curl -o ercc.txt https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt