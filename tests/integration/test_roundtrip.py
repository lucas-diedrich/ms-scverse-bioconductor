"""Round-trip integration tests between mulink/mudata and QFeatures.

Each test drives one full round trip and compares the two ends *in the language
the trip started in*, because "lossless" means lossless to that ecosystem's
object. Comparing the two intermediate .h5mu files instead would only test one
writer against itself.

    Python -> R -> Python    mudata writes, R reads and writes, mudata reads,
                             `compare_mudata()` compares MuData to MuData.

    R -> Python -> R         R writes, mudata reads and writes, R reads,
                             `compare_qfeatures()` compares QFeatures to
                             QFeatures and reports back as JSON.

Every case declares the deviation kinds it tolerates, and the assertion is an
equality: an unexpected kind is a regression, and a tolerated kind that stopped
appearing means the profile claims a loss that has since been fixed. That is
what stops the known-loss list from rotting into a list of excuses.

A case that cannot complete the trip at all is marked `xfail(strict=True)`
rather than given a tolerance, since there is no comparison to tolerate
anything about. Those marks flip to a failure the moment the underlying bug is
fixed, for the same reason.
"""

from __future__ import annotations

import anndata
import mudata as md
import numpy as np
import pytest
from anndata import AnnData

from ms_scverse_bioconductor.conformance import Report, compare_mudata

pytestmark = pytest.mark.integration


# Deviations tolerated on a `Python -> R -> Python` trip, per fixture.
#
# `set_order`/`global_feature_names` on `unsorted_modality_names`: modality
# order is lost in both directions. MuData (R) writes and reads it as the
# `order` attribute of the `/mod` group (`write_h5mu.R:180`, `utils.R:99`) while
# Python's mudata uses `mod-order`, so neither reads the other's and both fall
# back to HDF5's alphabetical name order. Because the global feature axis is the
# per-modality blocks concatenated in modality order, reordering the modalities
# also reorders `var_names`. The feature graph itself survives, which is why
# `feature_graph` is absent here: it is compared as a name-keyed edge set.
PYTHON_ORIGIN_PROFILE: dict[str, set[str]] = {
    "simulated": set(),
    "simulated_deep": set(),
    "annotated": set(),
    "missing_values": set(),
    "unsorted_modality_names": {"set_order", "global_feature_names"},
}

PYTHON_ORIGIN_XFAIL: dict[str, str] = {}


# Deviations tolerated on an `R -> Python -> R` trip, per fixture.
# `Rscript R/roundtrip-cli.R list-fixtures` prints the registry.
#
# `layer_names` (every fixture): MuData stores an assay's primary matrix as .X,
# which carries no name, so `assayNames()` comes back as `''` instead of the
# original name -- `character(0)` becomes `''`, and `c('assay', 'aggcounts')`
# becomes `c('', 'aggcounts')`.
#
# `feature_annotation_types` (feat1, feat3, nullable_rowdata): `readH5MU()`
# returns character rowData columns as factors, and `cast_nullable_columns()`
# casts NA-bearing integer/logical columns to double on the way out. Both are
# listed as known losses in ROADMAP.md section 4.
#
# `feature_names` + `feature_annotation_columns` (feat2, feat3): these two
# fixtures have row names that collide across assays, so `feature_index()`
# prefixes them with the assay name and keeps the original in
# `rowData$mulink_feature_id`. `qfeatures_read_h5mu()` has no inverse for that,
# so the prefixed names and the helper column both survive into the recovered
# object.
#
# `link_topology` (feat3): a consequence of the same missing inverse. The edges
# are still there and still connect the same assays -- `link_parents` does not
# deviate -- but their endpoints carry the prefixed feature names.
#
# `link_fcol` (feat3): `writeQFeaturesH5MU()` does not write
# `uns['mulink']['assays']` yet, so `AssayLink@fcol` has nowhere to go and the
# reader substitutes the .varp key. Documented in `io.R` and ROADMAP.md.
#
# `set_order` (feat3): the modality-order gap described above, seen from the R
# side. Only feat3 has more than one assay *and* a non-alphabetical order.
R_ORIGIN_PROFILE: dict[str, set[str]] = {
    "feat1": {"layer_names", "feature_annotation_types"},
    "feat2": {"layer_names", "feature_names", "feature_annotation_columns"},
    "feat3": {
        "set_order",
        "layer_names",
        "feature_names",
        "feature_annotation_columns",
        "feature_annotation_types",
        "link_topology",
        "link_fcol",
    },
    "ft_na": set(),
    "nullable_rowdata": {"layer_names", "feature_annotation_types"},
}

R_ORIGIN_XFAIL: dict[str, str] = {
    "ft_na": (
        "writeQFeaturesH5MU() casts NA-bearing integer/logical *annotation* columns to "
        "double but not the assay matrix itself, so an integer assay carrying NA is "
        "written through MuData's nullable-integer encoding, whose int8 mask anndata "
        "rejects with 'mask should be boolean numpy array'. The file is unreadable from "
        "Python, so the round trip never reaches the comparison."
    ),
}


def _cases(profile: dict[str, set[str]], xfails: dict[str, str]) -> list:
    """Parametrisation for one direction, marking the cases that cannot run."""
    return [
        pytest.param(
            name,
            marks=pytest.mark.xfail(strict=True, reason=xfails[name]) if name in xfails else (),
        )
        for name in sorted(profile)
    ]


def _assert_profile(report: Report, tolerated: set[str], label: str) -> None:
    """Assert that the observed deviations are exactly the tolerated ones."""
    observed = report.kinds()

    unexpected = observed - tolerated
    resolved = tolerated - observed

    assert not unexpected and not resolved, (
        f"{label}\n"
        f"  unexpected deviations: {sorted(unexpected) or 'none'}\n"
        f"  tolerated but no longer observed: {sorted(resolved) or 'none'}\n"
        f"  full report:\n{report.describe()}"
    )


@pytest.mark.parametrize("fixture_name", _cases(PYTHON_ORIGIN_PROFILE, PYTHON_ORIGIN_XFAIL))
def test_python_origin_roundtrip(fixture_name: str, python_fixtures, r_cli, tmp_path) -> None:
    original = python_fixtures[fixture_name]()

    written = tmp_path / "python-origin.h5mu"
    recovered = tmp_path / "python-origin.r.h5mu"

    original.write_h5mu(written)
    r_cli.run("read-write", "--in", str(written), "--out", str(recovered))

    # The file just written is the baseline, not the in-memory object:
    # write_h5mu() calls update(), which can pull modality columns into the
    # global frames, and that is a property of mudata rather than of the
    # conversion under test.
    report = compare_mudata(md.read_h5mu(written), md.read_h5mu(recovered))
    _assert_profile(report, PYTHON_ORIGIN_PROFILE[fixture_name], f"Python -> R -> Python on '{fixture_name}'")


@pytest.mark.parametrize("fixture_name", _cases(R_ORIGIN_PROFILE, R_ORIGIN_XFAIL))
def test_r_origin_roundtrip(fixture_name: str, r_cli, tmp_path) -> None:
    written = tmp_path / "r-origin.h5mu"
    recovered = tmp_path / "r-origin.py.h5mu"

    r_cli.run("fixture", "--name", fixture_name, "--out", str(written))

    # The Python leg: disk -> memory -> disk, in one process.
    md.read_h5mu(written).write_h5mu(recovered)

    report = Report.from_dict(r_cli.run_json("compare", "--fixture", fixture_name, "--in", str(recovered)))
    _assert_profile(report, R_ORIGIN_PROFILE[fixture_name], f"R -> Python -> R on '{fixture_name}'")


@pytest.mark.xfail(
    strict=True,
    reason="MuData 1.14.0 reads /var/_index as a dataset, so it cannot read the "
    "nullable-string-array group anndata >= 0.13 writes for a pandas 3 str index. "
    "The suite works around this via the `legacy_string_encoding` fixture; when "
    "this test starts passing, drop that fixture.",
)
def test_r_reads_default_string_encoding(r_cli, tmp_path) -> None:
    """Whether the R leg can read an .h5mu written with anndata's own defaults.

    This is the one place the string-encoding constraint is stated as a fact
    about the two stacks rather than worked around. It is `strict` on purpose:
    once MuData reads the newer encoding the test passes, the suite fails, and
    the workaround gets removed instead of quietly outliving its reason.
    """
    written = tmp_path / "default-encoding.h5mu"

    with anndata.settings.override(allow_write_nullable_strings=None):
        md.MuData({"mod0": AnnData(np.zeros((2, 2), dtype=float))}).write_h5mu(written)

    r_cli.run("read-write", "--in", str(written), "--out", str(tmp_path / "out.h5mu"))


def test_r_comparator_selftest(r_cli) -> None:
    """The R comparator detects each deviation kind it can report.

    The counterpart of `tests/test_conformance.py` for the R half. Without it,
    `R_ORIGIN_PROFILE` could pass vacuously: a comparator that reported nothing
    would satisfy every empty tolerance in it.
    """
    report = r_cli.run_json("selftest")

    failed = [case for case in report["cases"] if not case["ok"]]
    assert report["ok"], "\n".join(
        f"  {case['name']}: expected {case['expected']}, observed {case['observed']}" for case in failed
    )


def test_python_fixture_registry_matches_profile(python_fixtures) -> None:
    """The Python fixture registry and the tolerated-deviation profile agree.

    Without this, adding a fixture to ``conftest.PYTHON_FIXTURES`` would
    silently go untested, since the parametrisation is driven by the profile.
    """
    assert set(python_fixtures) == set(PYTHON_ORIGIN_PROFILE), (
        f"only in the registry: {sorted(set(python_fixtures) - set(PYTHON_ORIGIN_PROFILE))}; "
        f"only in the profile: {sorted(set(PYTHON_ORIGIN_PROFILE) - set(python_fixtures))}"
    )


def test_r_fixture_registry_matches_profile(r_cli) -> None:
    """The R fixture registry covers every fixture the profile names."""
    available = set(r_cli.run_json("list-fixtures")["fixtures"])
    assert set(R_ORIGIN_PROFILE) <= available, (
        f"profile names the fixtures {sorted(set(R_ORIGIN_PROFILE) - available)}, which R does not provide"
    )
