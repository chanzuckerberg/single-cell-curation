import logging
import sys

import click

logger = logging.getLogger("cellxgene_schema")


@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate the cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
@click.option("-v", "--verbose", help="When present will set logging level to debug", is_flag=True)
def schema_cli(verbose):
    logging.basicConfig(level=logging.ERROR)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)


@schema_cli.command(
    name="validate",
    short_help="Check that an h5ad follows the cellxgene data integration schema.",
    help="Check that an h5ad follows the cellxgene data integration schema. If validation fails this command will "
    "return an exit status of 1 otherwise 0. When the '--add-labels <FILE>' tag is present, the command will add "
    "ontology/gene labels based on IDs and write them to a new h5ad.",
)
@click.argument("h5ad_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-a",
    "--add-labels",
    "add_labels_file",
    help="When present it will add labels to genes and ontologies based on IDs",
    required=False,
    default=None,
    type=click.Path(exists=False, dir_okay=False, writable=True),
)
@click.option("-i", "--ignore-labels", help="Ignore ontology labels when validating", is_flag=True)
def schema_validate(h5ad_file, add_labels_file, ignore_labels):
    # Imports are very slow so we defer loading until Click arg validation has passed
    logger.info("Loading dependencies")
    try:
        import anndata  # noqa: F401
    except ImportError:
        raise click.ClickException("[cellxgene] cellxgene-schema requires anndata") from None

    logger.info("Loading validator modules")
    from .validate import validate

    is_valid, _, _ = validate(h5ad_file, add_labels_file, ignore_labels=ignore_labels)
    if is_valid:
        sys.exit(0)
    else:
        sys.exit(1)


@schema_cli.command(
    name="add-labels",
    short_help="Create a copy of an h5ad with portal-added labels",
    help="Create a copy of an h5ad with portal-added labels. The labels are added based on the IDs in the "
    "provided file.",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
def add_labels(input_file, output_file):
    from utils import read_h5ad

    from .write_labels import AnnDataLabelAppender

    logger.info(f"Loading h5ad from {input_file}")
    adata = read_h5ad(input_file)
    anndata_label_adder = AnnDataLabelAppender(adata)
    logger.info("Adding labels")
    if anndata_label_adder.write_labels(output_file):
        sys.exit(0)
    else:
        sys.exit(1)


@schema_cli.command(
    name="remove-labels",
    short_help="Create a copy of an h5ad without portal-added labels",
    help="Create a copy of an h5ad without portal-added labels.",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
def remove_labels(input_file, output_file):
    from .remove_labels import AnnDataLabelRemover

    logger.info("Loading dependencies")
    try:
        import anndata  # noqa: F401
    except ImportError:
        raise click.ClickException("[cellxgene] cellxgene-schema requires anndata") from None

    logger.info(f"Loading h5ad from {input_file}")
    adata = anndata.read_h5ad(input_file)
    anndata_label_remover = AnnDataLabelRemover(adata)
    if not anndata_label_remover.schema_def:
        return
    logger.info("Removing labels")
    anndata_label_remover.remove_labels()
    logger.info(f"Labels have been removed. Writing to {output_file}")
    anndata_label_remover.adata.write(output_file, compression="gzip")


@schema_cli.command(
    name="migrate",
    short_help="Convert an h5ad to the latest schema version.",
    help="Convert an h5ad from the previous to latest minor schema version. No validation will be "
    "performed on either the input or the output file.",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
@click.option("--collection_id", default=None, type=str, help="Collection ID, if migrating already published dataset")
@click.option("--dataset_id", default=None, type=str, help="Dataset ID, if migrating already published dataset")
def migrate(input_file, output_file, collection_id, dataset_id):
    from .migrate import migrate

    migrate(input_file, output_file, collection_id, dataset_id)


@click.command(
    name="map-species",
    short_help="Annotate non-human, non-mouse anndata with CL and UBERON equivalent terms",
    help="Annotate non-human, non-mouse anndata with CL and UBERON equivalent terms, based on values in"
    "organism-specific columns (e.g. organism_cell_type_ontology_term_id and organism_tissue_ontology_term_id)",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
def map_species(input_file, output_file):
    from .map_species import map_species

    map_species(input_file, output_file)


schema_cli.add_command(schema_validate)
schema_cli.add_command(migrate)
schema_cli.add_command(remove_labels)
schema_cli.add_command(map_species)

if __name__ == "__main__":
    schema_cli()
