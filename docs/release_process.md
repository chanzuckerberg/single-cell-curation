# cellxgene-schema Release Process

_This document defines the release process of the pip installable `cellxgene-schema` tool to PyPI._

The goal of the release process is to publish an installable package to PyPi, with a matching tagged release on github.

## Overview

The release process should result in the following side effects:

* Version number bump, using semantic versioning
* Tagged Github release
* Publication to PyPi

Note all release tags pushed to GitHub MUST follow semantic versioning.

## How to release the next version of `cellxgene-schema`

Steps must be run from the project directory and in a virtual env with all the dependencies defined in `requirements.txt` installed.

1. Decide on whether the changes going into the next release warrant a patch, minor, or major version bump using [semantic versioning](https://semver.org/) and determine the version number based on the current latest version (i.e. if the current version is 2.0.3 and you'd like to do a minor version bump, then the next version will be 2.1.0).

1. Create a release branch from latest main (e.g. `git checkout -b release-version-x.y.z`).

1. In `setup.py` change the version number to be the version number determined in step 1.

1. Update the [CHANGELOG](https://github.com/chanzuckerberg/single-cell-curation/blob/main/cellxgene_schema_cli/CHANGELOG.md) with the latest changes going into the new release. 

1. Commit the changes, push the new branch to origin, and open a PR to be reviewed. Do not merge in the PR yet.

1. Release the new version to Test PyPI. This can be done by running `make test-release`. Test the release by installing `cellxgene-schema` using pip `pip install -i https://test.pypi.org/simple/ cellxgene-schema` in a fresh virtual environment.

1. If you find errors with the release candidate, fix them in `main`, delete the current release branch, and start over from step 1 with the next version number. Please also update the CHANGELOG with a note that the version number has been burned. This is necessary because once a version goes out on Test PyPI (or PyPI), that version number cannot be used again (i.e. there are no updates supported). For example, if you were trying to release a minor bump, say 2.1.0, and you found an error, then the next version will be 2.2.0 and 2.1.0 will be a "burn".

1. If everything looks good, get the PR that you opened reviewed and **merged**. 

1. Release the version to (prod) PyPI by running `make release`.

1. Test the installation in a fresh virtual environment by running `pip install --no-cache-dir cellxgene-schema`. If you find errors, go back to step 1.

1. Create Github release using the version number and release notes ([instructions](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)):
    1. Draft new release. The release description should be the same as the CHANGELOG from Step 4. The release title should be "Release {version_num}".
    1. Select `main` as release branch (ensure you merged the release PR).
    1. Create new tag with same name as the release (`vx.y.z`).

## Gut-check checklist

High level checklist to ensure that you completed the release successfully.

- [ ] Do both the Test PyPI and PyPI latest versions match and are they the version you intended to release?
- [ ] Did you update the CHANGELOG with notes about what went into the new release?
- [ ] Did you create and merge in a PR to bump the version in `setup.py`?
- [ ] Did you create a Github release for the new version?
- [ ] Did you test the new release to make sure it works?

## Additional Resources

- [How to create a Python virtual environment](https://docs.python.org/3/library/venv.html#creating-virtual-environments)
- If you need access to the Test PyPI or PyPI projects, please ping the [#single-cell-eng](https://czi-sci.slack.com/archives/C023Q1APASK) Slack channel. 