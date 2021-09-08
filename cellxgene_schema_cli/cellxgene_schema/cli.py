import click
import sys
from cellxgene_schema import validate


@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate the cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
def schema_cli():
    try:
        import anndata  # noqa: F401
    except ImportError:
        raise click.ClickException("[cellxgene] cellxgene-schema requires anndata")


@click.command(
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
def schema_validate(h5ad_file, add_labels_file):
    if validate.validate(h5ad_file, add_labels_file):
        sys.exit(0)
    else:
        sys.exit(1)


schema_cli.add_command(schema_validate)

if __name__ == "__main__":
    schema_cli()
