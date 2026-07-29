"""Dump a mulink object to machine readable csv/tsv files"""

from pathlib import Path
from typing import Annotated

import networkx as nx
import pandas as pd
import typer
from mudata import MuData
from mulink.simulate import hierarchical_mudata

# Suffix follows the separator so that the files stay loadable by convention alone.
_SEPARATOR_EXTENSION = {"\t": "tsv", ",": "csv"}


def _extension(sep: str) -> str:
    return _SEPARATOR_EXTENSION.get(sep, "txt")


def _feature_edgelist(mdata: MuData, varp_key: str) -> pd.DataFrame:
    """Convert the feature mapping of ``mdata`` into a directed edge list.

    Parameters
    ----------
    mdata
        Object whose ``.varp[varp_key]`` holds the feature-mapping adjacency matrix.
    varp_key
        Key of the adjacency matrix in ``mdata.varp``.

    Returns
    -------
    Edge list with the columns ``source``, ``target``, ``weight``, ``source_mod``
    and ``target_mod``.
    """
    graph = nx.from_pandas_adjacency(
        pd.DataFrame.sparse.from_spmatrix(mdata.varp[varp_key], index=mdata.var_names, columns=mdata.var_names),
        create_using=nx.DiGraph,
    )

    edgelist = nx.to_pandas_edgelist(graph)

    # TODO: This relies on a convention in mulink
    # clarify between which layers the edges are drawn
    edgelist["source_mod"] = edgelist["source"].str.extract(r"(mod[0-9]+)")
    edgelist["target_mod"] = edgelist["target"].str.extract(r"(mod[0-9]+)")

    return edgelist


def export(
    n_mod: Annotated[int, typer.Option(help="Number of feature-levels in the object.")],
    n_obs: Annotated[int, typer.Option(help="Number of observations in the object.")] = 5,
    n_vertices: Annotated[int, typer.Option(help="Number of vertices in the level with the lowest cardinality.")] = 2,
    min_edges: Annotated[int, typer.Option(help="Minimum number of vertices between adjacent levels.")] = 2,
    extra_edge_probability: Annotated[
        float | None,
        typer.Option(
            help="Probability of adding an additional edge between vertices of adjacent levels. "
            "If not given, the feature relationship between levels is a tree, i.e. a feature of "
            "level n+1 maps to exactly one feature of level n."
        ),
    ] = 0.2,
    extra_edge_levels: Annotated[
        list[int] | None,
        typer.Option(
            help="Constrain the addition of extra edges to these levels (starting at 0). "
            "Repeat the option to pass several levels. Must not contain the highest level."
        ),
    ] = None,
    transitive_closure: Annotated[
        bool,
        typer.Option(
            help="Connect indirectly linked features with an explicit edge (A --> B --> C also yields A --> C)."
        ),
    ] = True,
    varp_key: Annotated[
        str, typer.Option(help="Key of the feature mapping in the `.varp` attribute of the simulated object.")
    ] = "feature_mapping",
    random_state: Annotated[int, typer.Option(help="Random state for the simulation.")] = 42,
    name: Annotated[str, typer.Option(help="Prefix of all written files.")] = "mudata",
    output_dir: Annotated[
        Path, typer.Option("--output-dir", "-o", help="Directory to write to. Created if it does not exist.")
    ] = Path(),
    sep: Annotated[
        str, typer.Option(help="Field separator of the written files. Determines the file extension.")
    ] = "\t",
    save_mudata: Annotated[
        bool, typer.Option("--save_mudata", "-save_mudata", help="Whether to save the MuData object as well")
    ] = True,
) -> None:
    """Simulate a hierarchical mulink object and dump it to machine readable text files.

    One file per omics layer holds the observation-by-feature matrix, one further
    file holds the feature mapping between the layers as directed edge list:

    ```
    <output-dir>/<name>.mod_<mod_name>.<ext>
    <output-dir>/<name>.edgelist.<ext>
    ```

    The extension is `tsv` for a tab separator, `csv` for a comma and `txt` otherwise.
    """
    mdata = hierarchical_mudata(
        n_mod=n_mod,
        n_obs=n_obs,
        n_vertices=n_vertices,
        min_edges=min_edges,
        extra_edge_probability=extra_edge_probability,
        extra_edge_levels=extra_edge_levels or None,
        transitive_closure=transitive_closure,
        varp_key=varp_key,
        random_state=random_state,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    extension = _extension(sep)

    for mod_name in mdata.mod:
        path = output_dir / f"{name}.mod_{mod_name}.{extension}"
        mdata[mod_name].to_df().to_csv(path, sep=sep, index_label="obs")
        typer.echo(f"Wrote {path}")

    path = output_dir / f"{name}.edgelist.{extension}"
    _feature_edgelist(mdata, varp_key=varp_key).to_csv(path, sep=sep, index=False)
    typer.echo(f"Wrote {path}")

    if save_mudata:
        path = output_dir / f"{name}.h5mu"
        mdata.write_h5mu(path)
        typer.echo(f"Wrote {path}")
