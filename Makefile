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

# RELEASE HELPERS
# For displaying the current version
current_version := $(shell awk '/current_version =/ { print $$3 }' .bumpversion.cfg)
pydist: clean
	cd cellxgene_schema_cli && python setup.py sdist bdist_wheel

# Show the current version
show-current-version:
	@echo $(current_version)

# Set PART={major,minor,patch}' as param to make bump.
# This will create a new release candidate. (i.e. 0.16.1 -> 0.16.2-rc.0 for a patch bump)
bump-version:
	bumpversion --config-file .bumpversion.cfg $(PART)

# This command is used to increment the rc number in case an initial candidate is erronous.
bump-release-candidate:
	bumpversion --config-file .bumpversion.cfg prerelversion --allow-dirty

# This command is used to remove the `rc` tag from the release candidate after it has been successfully tested in Test PyPI.
bump-release:
	bumpversion --config-file .bumpversion.cfg prerel --allow-dirty

# Create new version to commit to main
# Set PART=[major, minor, patch]
create-release-candidate: bump-version
	@echo "Version bumped part:$(PART) and client built. Ready to commit and push"

# Bump the release candidate version if needed (i.e. the previous release candidate had errors).
recreate-release-candidate: bump-release-candidate
	@echo "Version bumped part:$(PART) and client built. Ready to commit and push"

# Build dist and release the release candidate to Test PyPI
release-candidate-to-test-pypi: pydist
	pip install twine
	python -m twine upload --repository testpypi cellxgene_schema_cli/dist/*
	@echo "Release candidate dist built and uploaded to test.pypi.org"

# Build final dist (gets rid of the rc tag) and release final candidate to TestPyPI
release-final-to-test-pypi: bump-release pydist
	pip install twine
	python -m twine upload --repository testpypi cellxgene_schema_cli/dist/*
	@echo "Final release dist built and uploaded to test.pypi.org"
	@echo "Please test the release on Test PyPI"

release-final: pydist
	pip install twine
	python -m twine upload cellxgene_schema_cli/dist/*
	@echo "Release uploaded to pypi.org"
	@echo "Please test the release on PyPI"
