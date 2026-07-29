"""Command line interface of ms-scverse-bioconductor"""

import typer

from .export import export

# Markdown mode keeps the file-name schema in the help text un-rewrapped.
app = typer.Typer(add_completion=False, rich_markup_mode="markdown", no_args_is_help=True)

app.command()(export)


# A callback is what keeps typer from collapsing a single-command app into a
# top-level command, i.e. what makes `export` an actual subcommand.
@app.callback()
def main() -> None:
    """Helper scripts for the joint MS-proteomics analysis with Python and R packages."""


if __name__ == "__main__":
    app()
