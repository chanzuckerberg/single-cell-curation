.PHONY: update-references
update-references: download-ontologies gene-processing clean

.PHONY: download-ontologies
download-ontologies:
	python3 ./scripts/ontology_processing.py

.PHONY: gene-processing
gene-processing:
	python3 ./scripts/gene_processing.py

.PHONY: genes-update
genes-update:
	python3 . /scripts/update_genes.py

.PHONY: clean
clean:
	rm -rf temp
