# cellxgene-schema Release Process

_This document defines the release process of the pip installable `cellxgene-schema` tool to PyPI._

The goal of the release process is to publish an installable package to PyPI, with a matching tagged release on Github.

## Overview

The release process should result in the following side effects:

* Version number bump, using [semantic versioning](https://semver.org/)
* Tagged Github release
* Publication to PyPI

Note all release tags pushed to GitHub MUST follow semantic versioning.

## Prerequisites

Before attempting releases of `cellxgene-schema` to Test PyPI and PyPI, create an account on each:
* https://pypi.org/account/register/
* https://test.pypi.org/account/register/

and request to be granted the `Owner` role in the [#single-cell-eng](https://czi-sci.slack.com/archives/C023Q1APASK) Slack channel.
Once you have access, **please setup 2FA**! 

## How to release the next version of `cellxgene-schema`

Steps must be run from the project directory and in a virtual env with all the dependencies defined in `requirements.txt` installed.

1. Decide on whether the changes going into the next release warrant a patch, minor, or major version bump using [semantic versioning](https://semver.org/) and determine the version number based on the current latest version (i.e. if the current version is 2.0.4 and you'd like to do a minor version bump, then the next version will be 2.1.0).

1. Run `make create-release-candidate PART={major|minor|patch}`, where you set `PART` to be one of `major` or `minor` or `patch` based on Step 1. This step will automatically bump the version of the package to the appropriate next version in `setup.py` and commit the change to the branch. Note that this will create a release candidate version, so if the current version is 2.0.4 and you want to do a minor release, then the version will be `2.1.4-rc.0`. The `rc` suffix will be removed in a later step. This step will also automatically create and checkout a new release Git branch.

1. Update the [CHANGELOG](https://github.com/chanzuckerberg/single-cell-curation/blob/main/cellxgene_schema_cli/CHANGELOG.md) with the latest changes going into the new release. 

1. Commit the changes, push the new branch to origin, and open a PR to be reviewed. **Do not merge in the PR yet.**

1. Run `make release-candidate-to-test-pypi` to release the newly created release candidate to Test PyPI.

1. Using a fresh virtual environment, test the release candidate by installing `cellxgene-schema` using `python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ cellxgene-schema==2.1.0rc0`. (`--extra-index-url` is required to install dependencies that may not exist on Test PyPI)

1. If you detect errors in the release candidate, fix them in `main` and then rebase the release branch. Then run `make recreate-release-candidate` to bump up the release candidate version (i.e. 2.1.0-rc.0 will get bumped up to 2.1.0-rc.1). Run `make release-candidate-to-test-pypi` again to release this candidate to Test PyPI for testing.

1. If everything looks good, release the final candidate to Test PyPI by running `make release-final-to-test-pypi`. This will remove the `rc` tag from the version and upload the distribution to Test PyPI in one go. Test one last time to make sure everything is OK.

1. Ensure that all commits are pushed to your branch (especially make sure that all the version changes are pushed) and get the PR reviewed and **merged**.

1. Release to (prod) PyPI by running `make release-to-pypi`.

1. Test the installation in a fresh virtual environment by running `pip install --no-cache-dir cellxgene-schema`. If you find errors, go back to step 1.

1. Create a Github release ([instructions](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)):
    1. [Draft a new release](https://github.com/chanzuckerberg/single-cell-curation/releases/new).
      1. Create a new tag with same name as the release version (`x.y.z`).
      1. Select `main` as target branch (Ensure that you have merged the release PR).
      1. The release title should be `Release x.y.z`. The release description should be the same as the CHANGELOG from Step 4.

## Gut-check checklist

High level checklist to ensure that you completed the release successfully.

- [ ] Do both the Test PyPI and PyPI latest versions match and are they the version you intended to release?
- [ ] Did you update the CHANGELOG with notes about what went into the new release?
- [ ] Did you create and merge in a PR with the version in `setup.py` changed to the latest?
- [ ] Did you create a Github release for the new version?
- [ ] Did you test the new release to make sure it works?
- [ ] For any internal uses of `cellxgene-schema` have you upgraded the version number or let the appropriate stakeholders know that a new version exists?

## Additional Resources

- [How to create a Python virtual environment](https://docs.python.org/3/library/venv.html#creating-virtual-environments)