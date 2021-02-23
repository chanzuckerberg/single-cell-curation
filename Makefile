.PHONY: fmt
fmt:
	cd cellxgene_schema_cli/ && black cellxgene_schema/ tests/

lint:
	cd cellxgene_schema_cli/ && flake8 cellxgene_schema/ tests/

install:
	cd cellxgene_schema_cli/ && pip install -r requirements.txt && pip install .

unit-test:
	cd cellxgene_schema_cli && pytest

clean:
	rm -rf cellxgene_schema_cli/build cellxgene_schema_cli/dist cellxgene_schema_cli/cellxgene_schema.egg-info

pydist: clean
	cd cellxgene_schema_cli && python setup.py sdist bdist_wheel

test-release: pydist
	pip install twine
	python -m twine upload --repository testpypi cellxgene_schema_cli/dist/*

release: pydist
	pip install twine
	python -m twine upload cellxgene_schema_cli/dist/*
