.PHONY: update-references
update-references: download-ontologies gene-processing clean

.PHONY: download-ontologies
download-ontologies:
	python3 ./scripts/ontology_processing.py

.PHONY: gene-processing
gene-processing: download-gtf-10x-human download-gtf-10x-mouse download-gtf-covid19 download-gtf-ercc
	python3 ./scripts/gene_processing.py

.PHONY: download-gtf-10x-human
download-gtf-10x-human:
	mkdir -p temp
	echo Downloading human GTF
	curl -o ./temp/homo_sapiens.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

.PHONY: download-gtf-10x-mouse
download-gtf-10x-mouse:
	mkdir -p temp
	echo Downloading human GTF
	curl -o ./temp/mus_musculus.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

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
	rm -rf temp
