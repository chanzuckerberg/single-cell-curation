.PHONY: fmt
fmt:
	pre-commit run -a

lint:
	cd cellxgene_schema_cli/ && flake8 cellxgene_schema/ tests/

install:
	cd cellxgene_schema_cli/ && pip install -r requirements.txt && pip install .

unit-test:
	cd cellxgene_schema_cli && coverage run --source=cellxgene_schema -m pytest --log-level=INFO ./tests  \
		&& coverage xml

ontology-dry-run-tests:
	cd schema_bump_dry_run_scripts/ontologies/tests && pytest

clean:
	rm -rf cellxgene_schema_cli/build cellxgene_schema_cli/dist cellxgene_schema_cli/cellxgene_schema.egg-info

### RELEASE HELPERS ###

pydist: clean
	cd cellxgene_schema_cli && python setup.py sdist bdist_wheel

# For displaying the current version
current_version = $$(awk '/current_version =/ { print $$3 }' .bumpversion.cfg)

# Show the current version
show-current-version:
	@echo $(current_version)

# Create a release candidate version from a previously released version.  Commit then version change to a new release branch.
# Set PART={major,minor,patch} as param to make the version bump happen automatically.
create-release-candidate:
	bumpversion --config-file .bumpversion.cfg $(PART) --allow-dirty
	git checkout -b release-version-$$(awk '/current_version =/ { split($$3,v,"-"); print v[1] }' .bumpversion.cfg)
	@echo "Version bumped part:$(PART) to "$(current_version)". Ready to commit and push"

# Create another release candidate version from a previously created release candidate version, in case the previous release candidate had errors. Commit the version change to the current branch.
recreate-release-candidate:
	bumpversion --config-file .bumpversion.cfg prerelversion --allow-dirty
	@echo "Version bumped part:rc to "$(current_version)". Ready to commit and push"

# Build dist and release the release candidate to Test PyPI
release-candidate-to-test-pypi: pydist
	pip install twine
	python -m twine upload --repository testpypi cellxgene_schema_cli/dist/*
	@echo "Release candidate dist built and uploaded to test.pypi.org"

# Build final dist (gets rid of the rc tag) and release final candidate to TestPyPI
release-final-to-test-pypi:
	bumpversion --config-file .bumpversion.cfg prerel --allow-dirty --tag
	$(MAKE) pydist
	python -m twine upload --repository testpypi cellxgene_schema_cli/dist/*
	@echo "Final release dist built for "$(current_version)" and uploaded to test.pypi.org"
	@echo "Please test the release on Test PyPI"

release-to-pypi: pydist
	pip install twine
	python -m twine upload cellxgene_schema_cli/dist/*
	@echo "Release uploaded to pypi.org"
	@echo "Please test the release on PyPI"
