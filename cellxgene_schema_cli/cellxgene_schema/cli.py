import click
import sys


@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate the cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
def schema_cli():
    pass


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
@click.option(
    "-i",
    "--ignore-labels",
    help="Ignore ontology labels when validating",
    is_flag=True
)
@click.option(
    "-v",
    "--verbose",
    help="When present will set logging level to debug",
    is_flag=True
)
def schema_validate(h5ad_file, add_labels_file, ignore_labels, verbose):
    # Imports are very slow so we defer loading until Click arg validation has passed

    print("Loading dependencies")
    try:
        import anndata  # noqa: F401
    except ImportError:
        raise click.ClickException("[cellxgene] cellxgene-schema requires anndata")

    print("Loading validator modules")
    from .validate import validate

    is_valid, _, _ = validate(h5ad_file, add_labels_file, ignore_labels=ignore_labels, verbose=verbose)
    if is_valid:
        sys.exit(0)
    else:
        sys.exit(1)

@click.command(
    name="convert",
    short_help="Convert an h5ad from version 2.0.0 to version 3.0.0",
    help="Convert an h5ad from version 2.0.0 to version 3.0.0. No validation will be performed on either"
    "the input or the output file.",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
def convert(input_file, output_file):
    from .convert import convert
    convert(input_file, output_file)

schema_cli.add_command(schema_validate)
schema_cli.add_command(convert)

if __name__ == "__main__":
    schema_cli()
