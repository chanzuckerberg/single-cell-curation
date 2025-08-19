from setuptools import find_packages, setup

with open("requirements.txt") as fh:
    requirements = fh.read().splitlines()

setup(
    name="cellxgene-schema",
    version="6.0.1",
    url="https://github.com/chanzuckerberg/single-cell-curation",
    license="MIT",
    author="Chan Zuckerberg Initiative",
    author_email="cellxgene@chanzuckerberg.com",
    description="Tool for applying and validating cellxgene integration schema to single cell datasets",
    long_description="Tool for applying and validating cellxgene integration schema to single cell datasets",
    install_requires=requirements,
    python_requires=">=3.10",
    packages=find_packages(exclude=["tests", "scripts"]),
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={"console_scripts": ["cellxgene-schema = cellxgene_schema.cli:schema_cli"]},
)
