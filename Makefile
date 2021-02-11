.PHONY: fmt
fmt:
	black cellxgene_schema/ tests/

lint:
	flake8 cellxgene_schema/ tests/

install:
	pip install -r requirements.txt
	pip install .

unit-test:
	pytest
