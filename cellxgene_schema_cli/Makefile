.PHONY: update-references
update-references: download-ontologies gene-processing clean

.PHONY: download-ontologies
download-ontologies:
	python3 ./scripts/ontology_processing.py

.PHONY: gene-processing
gene-processing: download-gtf-human download-gtf-mouse download-gtf-covid19 download-gtf-ercc
	python3 ./scripts/gene_processing.py

.PHONY: download-gtf-human
download-gtf-human:
	mkdir -p temp
	echo Downloading human GTF
	curl -o ./temp/homo_sapiens.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz

.PHONY: download-gtf-mouse
download-gtf-mouse:
	mkdir -p temp
	echo Downloading mouse GTF
	curl -o ./temp/mus_musculus.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz

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
