import click

from cellxgene_schema import ontology

@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate the cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)

def schema_cli():
    try:
        import anndata
    except ImportError:
        raise click.ClickException("[cellxgene] cellxgene schema requires anndata")

@click.command(
    name="placeholder",
)
def placeholder():
    click.echo("placeholder")


schema_cli.add_command(placeholder)
