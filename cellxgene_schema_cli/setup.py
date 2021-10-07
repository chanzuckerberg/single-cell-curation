from setuptools import setup, find_packages

with open("requirements.txt") as fh:
    requirements = fh.read().splitlines()

setup(
    name="cellxgene-schema",
    version="2.0.2",
    url="https://github.com/chanzuckerberg/single-cell-curation",
    license="MIT",
    author="Chan Zuckerberg Initiative",
    author_email="cellxgene@chanzuckerberg.com",
    description="Tool for applying and validating cellxgene integration schema to single cell datasets",
    install_requires=requirements,
    packages=["cellxgene_schema"],
    package_dir={"cellxgene_schema": "cellxgene_schema"},
    package_data={"cellxgene_schema": ["ontology_files/*gz", "schema_definitions/*yaml"]},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={"console_scripts": ["cellxgene-schema = cellxgene_schema.cli:schema_cli"]}
)
