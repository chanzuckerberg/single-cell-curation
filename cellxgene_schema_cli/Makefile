.PHONY: update-gene-ontologies
update-gene-ontologies: gene-processing clean

.PHONY: gene-processing
gene-processing:
	python3 ./scripts/gene_processing.py

.PHONY: genes-update
genes-update:
	python3 ./scripts/update_genes.py

.PHONY: clean
clean:
	rm -rf temp
