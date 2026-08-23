"""Unit tests for the Python conformance comparator.

These need no R. They pin that the comparator actually detects a change, which
is the property the round-trip suite depends on: a comparator that reports
nothing would make every profile in `tests/integration` pass vacuously.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from mudata import MuData
from scipy.sparse import csr_matrix

from ms_scverse_bioconductor.conformance import Deviation, Report, compare_mudata

VARP_KEY = "feature_mapping"


def _modality(name: str, n_var: int, n_obs: int = 4) -> ad.AnnData:
    rng = np.random.default_rng(abs(hash(name)) % 2**32)
    return ad.AnnData(
        X=rng.random((n_obs, n_var)),
        var=pd.DataFrame({"level": [name] * n_var}, index=[f"{name}-{index}" for index in range(n_var)]),
        obs=pd.DataFrame(
            {"batch": [f"b{index % 2}" for index in range(n_obs)]}, index=[f"c{index}" for index in range(n_obs)]
        ),
    )


def _build(order: tuple[str, ...] = ("parent", "child")) -> MuData:
    """A two-modality object with a small parent -> child feature graph."""
    mods = {"parent": _modality("parent", 3), "child": _modality("child", 4)}
    mdata = MuData({name: mods[name] for name in order})

    positions = {name: index for index, name in enumerate(mdata.var_names)}
    edges = [("parent-0", "child-0"), ("parent-0", "child-1"), ("parent-1", "child-2"), ("parent-2", "child-3")]
    rows = [positions[source] for source, _ in edges]
    cols = [positions[target] for _, target in edges]

    mdata.varp[VARP_KEY] = csr_matrix((np.ones(len(edges)), (rows, cols)), shape=(mdata.n_vars, mdata.n_vars))
    return mdata


def test_identical_objects_report_no_deviation() -> None:
    report = compare_mudata(_build(), _build())

    assert report.equal
    assert report.kinds() == set()
    assert report.describe() == "no deviations"


def test_changed_value_is_reported() -> None:
    recovered = _build()
    recovered.mod["child"].X[0, 0] += 1.0

    assert compare_mudata(_build(), recovered).kinds() == {"values"}


def test_introduced_missing_value_is_reported() -> None:
    recovered = _build()
    recovered.mod["child"].X[1, 1] = np.nan

    assert compare_mudata(_build(), recovered).kinds() == {"values"}


def test_equal_missingness_pattern_is_not_reported() -> None:
    """NaN in the same place on both sides is agreement, not a deviation."""
    original, recovered = _build(), _build()
    original.mod["child"].X[1, 1] = np.nan
    recovered.mod["child"].X[1, 1] = np.nan

    assert compare_mudata(original, recovered).equal


def test_renamed_feature_is_reported() -> None:
    """A modality's own feature names are compared independently of the global axis.

    Renaming a modality's ``var_names`` without calling ``update()`` leaves the
    global ``var_names`` stale, so only the per-modality deviation is expected
    here; the global axis has its own kind.
    """
    recovered = _build()
    names = list(recovered.mod["child"].var_names)
    names[0] = "renamed"
    recovered.mod["child"].var_names = names

    assert compare_mudata(_build(), recovered).kinds() == {"feature_names"}


def test_reordered_modalities_are_reported_as_ordering() -> None:
    """A reordering must not masquerade as a corrupted feature graph.

    The graph is compared as a name-keyed edge set precisely so that moving the
    global feature axis reports as `set_order`/`global_feature_names` and leaves
    `feature_graph` alone.
    """
    report = compare_mudata(_build(("parent", "child")), _build(("child", "parent")))

    assert report.kinds() == {"set_order", "global_feature_names"}


def test_added_annotation_column_is_reported() -> None:
    recovered = _build()
    recovered.mod["child"].var["extra"] = 1.0

    assert compare_mudata(_build(), recovered).kinds() == {"feature_annotation_columns"}


def test_changed_annotation_dtype_is_reported_as_a_type_change() -> None:
    """Same values, different storage type: a type deviation and nothing else."""
    original, recovered = _build(), _build()
    original.mod["child"].var["count"] = np.arange(original.mod["child"].n_vars, dtype="int64")
    recovered.mod["child"].var["count"] = np.arange(recovered.mod["child"].n_vars, dtype="float64")

    assert compare_mudata(original, recovered).kinds() == {"feature_annotation_types"}


def test_dropped_graph_edge_is_reported() -> None:
    recovered = _build()
    graph = recovered.varp[VARP_KEY].tolil()
    graph[recovered.var_names.get_loc("parent-0"), recovered.var_names.get_loc("child-0")] = 0
    recovered.varp[VARP_KEY] = graph.tocsr()

    assert compare_mudata(_build(), recovered).kinds() == {"feature_graph"}


def test_missing_modality_is_reported() -> None:
    original = _build()
    recovered = MuData({"parent": original.mod["parent"].copy()})

    report = compare_mudata(original, recovered)
    assert "sets" in report.kinds()


def test_changed_shape_is_reported_without_amplification() -> None:
    """A shape mismatch is reported alone, since everything else is positional."""
    original = _build()
    recovered = MuData({"parent": original.mod["parent"].copy(), "child": _modality("child", 2)})

    report = compare_mudata(original, recovered)
    assert "shape" in report.kinds()
    assert {deviation.kind for deviation in report.deviations if deviation.scope == "child"} == {"shape"}


@pytest.mark.parametrize("payload", [{}, {"deviations": []}])
def test_report_from_empty_payload_is_equal(payload: dict) -> None:
    assert Report.from_dict(payload).equal


def test_report_from_dict_matches_the_r_comparator_shape() -> None:
    payload = {
        "equal": False,
        "deviations": [
            {"kind": "link_fcol", "scope": "peptides", "detail": "a: [Sequence]; b: [feature_mapping]"},
            {"kind": "set_order", "scope": "", "detail": "a: [x,y]; b: [y,x]"},
        ],
    }

    report = Report.from_dict(payload)

    assert not report.equal
    assert report.kinds() == {"link_fcol", "set_order"}
    assert report.deviations[0] == Deviation("link_fcol", "peptides", "a: [Sequence]; b: [feature_mapping]")
    assert "link_fcol: [peptides]" in report.describe()
