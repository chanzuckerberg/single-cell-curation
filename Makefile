.PHONY: fmt
fmt:
	cd cellxgene_schema_cli/ && black cellxgene_schema/ tests/

lint:
	cd cellxgene_schema_cli/ && flake8 cellxgene_schema/ tests/

install:
	cd cellxgene_schema_cli/ && pip install -r requirements.txt && pip install .

unit-test:
	cd cellxgene_schema_cli && pytest
