"""Fixtures for the mulink <-> QFeatures round-trip integration tests.

The two stacks run as separate processes and exchange .h5mu files only. This
module owns the R side of that contract: locating an R interpreter that can
actually import QFeatures and MuData, and skipping the suite rather than failing
it when none is available, so the unit suite stays green on a machine without R.

Set ``MSSB_RSCRIPT`` to point at a specific interpreter, for instance the one in
a mamba environment::

    MSSB_RSCRIPT=~/mamba/envs/qfeatures/bin/Rscript pytest tests/integration
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import anndata
import anndata as ad
import numpy as np
import pytest
from mudata import MuData
from mulink.simulate import hierarchical_mudata

VARP_KEY = "feature_mapping"

#: Repository root, from ``tests/integration/conftest.py``.
ROOT = Path(__file__).resolve().parents[2]

R_CLI = ROOT / "R" / "roundtrip-cli.R"


def pytest_configure(config: pytest.Config) -> None:
    """Register the markers this suite uses."""
    config.addinivalue_line("markers", "integration: requires both the R and the Python stack")
    config.addinivalue_line("markers", "slow: downloads or processes a realistic dataset")


@dataclass(frozen=True)
class RCli:
    """Thin wrapper around ``Rscript R/roundtrip-cli.R``."""

    executable: str
    script: Path

    def run(self, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
        """Invoke a subcommand, raising with R's stderr attached on failure."""
        completed = subprocess.run(
            [self.executable, str(self.script), *args],
            capture_output=True,
            text=True,
            cwd=cwd,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"`{' '.join(args)}` failed with exit code {completed.returncode}\n"
                f"--- stdout ---\n{completed.stdout}\n--- stderr ---\n{completed.stderr}"
            )
        return completed

    def run_json(self, *args: str, cwd: Path | None = None) -> dict:
        """Invoke a subcommand whose stdout is a JSON document."""
        return json.loads(self.run(*args, cwd=cwd).stdout)


@pytest.fixture(scope="session")
def r_cli() -> RCli:
    """An R interpreter that can run the round-trip CLI, or skip the test."""
    executable = os.environ.get("MSSB_RSCRIPT") or shutil.which("Rscript")
    if executable is None:
        pytest.skip("No Rscript on PATH; set MSSB_RSCRIPT to select an interpreter.")

    cli = RCli(executable=executable, script=R_CLI)

    # `check` reports the missing packages instead of dying inside a library()
    # call, which is the difference between an actionable skip and a traceback.
    try:
        completed = subprocess.run(
            [executable, str(R_CLI), "check"],
            capture_output=True,
            text=True,
        )
    except OSError as error:
        # Reached when MSSB_RSCRIPT names a path that is not executable. Skipping
        # keeps the suite non-blocking, and naming the path makes the
        # misconfiguration obvious rather than looking like a missing stack.
        pytest.skip(f"Cannot execute '{executable}': {error}")

    if completed.returncode != 0:
        try:
            missing = ", ".join(json.loads(completed.stdout).get("missing", []))
        except json.JSONDecodeError:
            missing = completed.stderr.strip() or "unknown"
        pytest.skip(f"{executable} cannot run the R leg; missing: {missing}")

    return cli


@pytest.fixture(scope="session", autouse=True)
def legacy_string_encoding():
    """Write string arrays in the encoding MuData (R) is able to read.

    pandas 3 gives an index the new ``str`` dtype, and anndata >= 0.13 stores
    such a column as ``nullable-string-array``: a *group* holding ``values`` and
    ``mask``. MuData 1.14.0 reads ``/var/_index`` as a *dataset*, so it fails
    with "The provided H5Identifier is not a dataset identifier" and the whole
    exchange stops at the index, before any of the conversion under test runs.

    ``allow_write_nullable_strings = False`` restores the ``string-array``
    encoding, which both sides read. The constraint is asserted separately by
    ``test_r_reads_default_string_encoding``, so it disappears from here the
    moment MuData learns the newer encoding rather than lingering unnoticed.
    """
    if not hasattr(anndata.settings, "allow_write_nullable_strings"):
        yield
        return

    previous = anndata.settings.allow_write_nullable_strings
    anndata.settings.allow_write_nullable_strings = False
    try:
        yield
    finally:
        anndata.settings.allow_write_nullable_strings = previous


def _annotate(adata: ad.AnnData, name: str, rng: np.random.Generator) -> ad.AnnData:
    """Attach the annotation types a converter has to carry across languages."""
    adata.var["level"] = name
    adata.var["score"] = rng.random(adata.n_vars)
    adata.obs["batch"] = [f"batch{index % 2}" for index in range(adata.n_obs)]
    adata.obs["depth"] = rng.random(adata.n_obs)
    return adata


def _rebuild(mdata: MuData, mods: dict[str, ad.AnnData]) -> MuData:
    """Rebuild a MuData from replacement modalities, keeping the feature axis.

    The global feature axis is the per-modality blocks concatenated in modality
    order, so as long as the insertion order and the per-modality feature names
    are preserved, the axis is too and the ``.varp`` matrix stays aligned. The
    assertion pins that, because ``.varp`` is meaningless if it silently stops
    matching ``var_names``.
    """
    graph = mdata.varp[VARP_KEY]
    rebuilt = MuData(mods)

    assert list(rebuilt.var_names) == list(mdata.var_names), "rebuilding reordered the global feature axis"
    rebuilt.varp[VARP_KEY] = graph
    return rebuilt


def _simulated() -> MuData:
    """Three feature levels, a tree-shaped graph plus a few extra edges."""
    return hierarchical_mudata(n_mod=3, n_obs=8, n_vertices=3, transitive_closure=False, random_state=0)


def _simulated_deep() -> MuData:
    """Four feature levels with denser fan-in between adjacent levels."""
    return hierarchical_mudata(
        n_mod=4,
        n_obs=6,
        n_vertices=2,
        extra_edge_probability=0.5,
        transitive_closure=False,
        random_state=1,
    )


def _annotated() -> MuData:
    """Three feature levels carrying `.var` and `.obs` annotation columns."""
    rng = np.random.default_rng(2)
    mdata = hierarchical_mudata(n_mod=3, n_obs=8, n_vertices=3, transitive_closure=False, random_state=2)
    return _rebuild(mdata, {name: _annotate(mdata.mod[name].copy(), name, rng) for name in mdata.mod})


def _missing_values() -> MuData:
    """Three feature levels with an NaN pattern in every `.X`."""
    rng = np.random.default_rng(3)
    mdata = hierarchical_mudata(n_mod=3, n_obs=8, n_vertices=3, transitive_closure=False, random_state=3)
    for name in mdata.mod:
        adata = mdata.mod[name]
        mask = rng.random(adata.shape) < 0.25
        adata.X = np.where(mask, np.nan, adata.X)
    return mdata


def _unsorted_modality_names() -> MuData:
    """Modality names whose insertion order is not their alphabetical order.

    ``precursors, peptides, proteins`` is the real feature hierarchy and sorts
    to ``peptides, precursors, proteins``, so this fixture is what makes a lost
    modality order observable from the Python side; every ``mod<n>`` fixture is
    already in alphabetical order and cannot show it.
    """
    mdata = hierarchical_mudata(n_mod=3, n_obs=8, n_vertices=3, transitive_closure=False, random_state=4)
    mapping = {"mod0": "precursors", "mod1": "peptides", "mod2": "proteins"}
    return _rebuild(mdata, {mapping[name]: mdata.mod[name].copy() for name in mdata.mod})


#: Python-origin fixtures, by name. Each is deterministic.
PYTHON_FIXTURES: dict[str, Callable[[], MuData]] = {
    "simulated": _simulated,
    "simulated_deep": _simulated_deep,
    "annotated": _annotated,
    "missing_values": _missing_values,
    "unsorted_modality_names": _unsorted_modality_names,
}


@pytest.fixture(scope="session")
def python_fixtures() -> dict[str, Callable[[], MuData]]:
    """The registry of Python-origin fixture builders."""
    return PYTHON_FIXTURES
