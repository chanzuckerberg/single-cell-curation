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


def _check_anndata_requires_fragment(h5ad_file):
    from .atac_seq import check_anndata_requires_fragment, report_errors

    try:
        fragment_required = check_anndata_requires_fragment(h5ad_file)
        if fragment_required:
            logger.info("Anndata requires an ATAC fragment file.")
        else:
            logger.info("Anndata does not require an ATAC fragment file.")
    except Exception as e:
        report_errors("Anndata does not support ATAC fragment files for the followings reasons:", [str(e)])
        sys.exit(1)


@schema_cli.command(
    name="process-fragment",
    short_help="Check that an ATAC SEQ fragment follows the cellxgene data integration schema.",
    help="Check that an ATAC SEQ fragment follows the cellxgene data integration schema. If validation fails this "
    "command will return an exit status of 1 otherwise 0. When the '--generate-index' tag is present, "
    "the command will generate a tabix compatible version of the fragment and tabix index. The generated "
    "fragment will have the file suffix .bgz and the index will have the file suffix .bgz.tbi.",
)
@click.argument("h5ad_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("fragment_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.option("-i", "--generate-index", help="Generate index for fragment", is_flag=True)
@click.option("-o", "--output-file", help="Output file for the processed fragment.", type=click.Path(exists=False))
def fragment_validate(h5ad_file, fragment_file, generate_index, output_file):
    from .atac_seq import process_fragment

    _check_anndata_requires_fragment(h5ad_file)

    if not process_fragment(fragment_file, h5ad_file, generate_index=generate_index, output_file=output_file):
        sys.exit(1)


# add a cli command to deduplicate an ATAC fragment file
@schema_cli.command(
    name="deduplicate-fragment",
    short_help="Deduplicate an ATAC SEQ fragment file.",
    help="Deduplicate an ATAC SEQ fragment file. If deduplication fails this command will return an exit status of 1 "
    "otherwise 0. The deduplicated fragment will have the file suffix .dedup.bgz and the index will have the file "
    "suffix .dedup.bgz.tbi.",
)
@click.argument("fragment_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.option("-o", "--output-file", help="Output file for the deduplicated fragment.", type=click.Path(exists=False))
@click.option("-m", "--memory", help="Memory limit as a percentage of total memory.", type=int, default=80)
def deduplicate_fragment(fragment_file, output_file, memory):
    from .atac_seq import deduplicate_fragment_rows

    try:
        deduplicate_fragment_rows(fragment_file, output_file_name=output_file, sort_memory_percent=memory)
    except Exception as e:
        logger.error(f"Failed to deduplicate fragment file: {e}")
        sys.exit(1)


@schema_cli.command(
    name="check-anndata-requires-fragment",
    short_help="Check if that the anndata provided supports an Atac seq fragment file.",
    help="Check that an ATAC SEQ anndata.obs['assay_ontology_term_id'] is all paired or unpaired assays. "
    "This determines if a fragment file is required, optional, or forbiden. "
    "If the anndata does not support a fragment, an error message will be printed and exit status will be 1.",
)
@click.argument("h5ad_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
def check_anndata_requires_fragment(h5ad_file):
    _check_anndata_requires_fragment(h5ad_file)


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


@schema_cli.command(
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


if __name__ == "__main__":
    schema_cli()
