.PHONY: update-references
update-references: download-ontologies gene_processing clean

.PHONY: download-ontologies
download-ontologies:
	python3 ./scripts/ontology_processing.py

.PHONY: gene-processing
gene-processing: download-gtf-10x-human download-gtf-10x-mouse download-gtf-covid19 download-gtf-ercc
	python3 ./scripts/gene_processing.py

.PHONY: download-gtf-10x-human
download-gtf-10x-human:
	$(eval url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz")
	$(eval base_dir="./temp/refdata-gex-GRCh38-2020-A")
	$(eval output_file="./temp/homo_sapiens.gtf.gz")
	$(MAKE) url=$(url) base_dir=$(base_dir) output_file=$(output_file) download-gtf-10x

.PHONY: download-gtf-10x-mouse
download-gtf-10x-mouse:
	$(eval url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz")
	$(eval base_dir="./temp/refdata-gex-mm10-2020-A")
	$(eval output_file="./temp/mus_musculus.gtf.gz")
	$(MAKE) url=$(url) base_dir=$(base_dir) output_file=$(output_file) download-gtf-10x

.PHONY: download-gtf-10x
download-gtf-10x:
	mkdir -p temp
	cd temp

	## Download data and copies data
	curl -o temp/10Xdata.tar.gz ${url}
	tar -xvzf temp/10Xdata.tar.gz -C temp

	## Rename base dir because it changes from human to mouse
	mv temp/${base_dir} temp/10Xdata
	cp temp/10Xdata/genes/genes.gtf temp/${output_file}
	rm -r temp/10Xdata

.PHONY: download-gtf-covid19
download-gtf-covid19:
	mkdir -p temp
	echo Downloading sars_cov_2 GTF
	curl -o temp/sars_cov_2.gtf.gz ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz

.PHONY: download-gtf-ercc
download-gtf-ercc:
	mkdir -p temp
	echo Downloading ERCC gene ids
	curl -o ./temp/ercc.txt https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt

.PHONY: clean
clean:
	rm -r temp