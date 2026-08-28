"""Conformance comparison of two MuData objects.

A comparison returns a list of deviations rather than a single boolean, because
"lossless" is a profile and not a property: some deviations are known and
accepted, and a test states which ones it tolerates. Reporting kinds separately
is what keeps that list honest- once a gap is closed the kind disappears and
the test that still tolerates it fails.

A round trip must be compared in the language it started in, since "lossless"
means lossless to *that* ecosystem's object. So this module compares the two
ends of a ``Python -> R -> Python`` trip, while ``R/compare.R`` compares the two
ends of an ``R -> Python -> R`` trip. Comparing the intermediate .h5mu files
instead would test one writer against itself.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
from scipy.sparse import issparse, spmatrix

#: Deviation kinds. Shared with the R comparator; only the subset that applies
#: to a given language is ever emitted.
#:
#: ``sets``                        set of modality/assay names differs
#: ``set_order``                   same names, different order
#: ``layer_names``                 ``.layers`` of a modality differs
#: ``shape``                       modality dimensions differ
#: ``feature_names``               per-modality ``var_names`` differ
#: ``obs_names``                   per-modality ``obs_names`` differ
#: ``values``                      ``.X`` values or missingness differ
#: ``feature_annotation_*``        per-modality ``.var`` columns/types/values
#: ``obs_annotation_*``            per-modality ``.obs`` columns/types/values
#: ``global_feature_names``        global ``var_names`` -- the ``.varp`` axis
#: ``global_obs_names``            global ``obs_names``
#: ``global_obs_annotation_*``     global ``.obs`` columns/types/values
#: ``feature_graph``               ``.varp`` edge set
#: ``feature_graph_keys``          set of ``.varp`` keys
DEVIATION_KINDS = (
    "sets",
    "set_order",
    "layer_names",
    "shape",
    "feature_names",
    "obs_names",
    "values",
    "feature_annotation_columns",
    "feature_annotation_types",
    "feature_annotation_values",
    "obs_annotation_columns",
    "obs_annotation_types",
    "obs_annotation_values",
    "global_feature_names",
    "global_obs_names",
    "global_obs_annotation_columns",
    "global_obs_annotation_types",
    "global_obs_annotation_values",
    "feature_graph",
    "feature_graph_keys",
)

_TOLERANCE = 1e-8


@dataclass(frozen=True)
class Deviation:
    """One way in which a round trip changed the object.

    Attributes
    ----------
    kind
        One of :data:`DEVIATION_KINDS`.
    scope
        The modality the deviation applies to, or ``""`` for the object as a
        whole.
    detail
        Human readable description, for the test failure message.
    """

    kind: str
    scope: str
    detail: str

    def __str__(self) -> str:
        scope = f"[{self.scope}] " if self.scope else ""
        return f"{self.kind}: {scope}{self.detail}"


@dataclass
class Report:
    """The outcome of one conformance comparison."""

    deviations: list[Deviation] = field(default_factory=list)

    @property
    def equal(self) -> bool:
        """Whether the two objects agreed on every compared property."""
        return not self.deviations

    def kinds(self) -> set[str]:
        """The distinct deviation kinds observed, for comparison to a profile."""
        return {deviation.kind for deviation in self.deviations}

    def describe(self) -> str:
        """Multi-line rendering of every deviation, for a failure message."""
        if self.equal:
            return "no deviations"
        return "\n".join(f"  - {deviation}" for deviation in self.deviations)

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> Report:
        """Rebuild a report from the JSON emitted by the R comparator."""
        return cls(
            [
                Deviation(record["kind"], record.get("scope", ""), record.get("detail", ""))
                for record in payload.get("deviations", [])
            ]
        )


def _dense(matrix: np.ndarray | spmatrix) -> np.ndarray:
    return np.asarray(matrix.todense() if issparse(matrix) else matrix)


def _equal_matrix(a: Any, b: Any) -> bool:
    """Compare two quantitative matrices, tolerating NA/NaN interchange.

    Missingness is compared as a pattern and values only where both are
    present. R's ``NA`` reaches Python as ``NaN``, so the two cannot be told
    apart across an .h5mu round trip; that is a documented known loss rather
    than something to report per matrix.
    """
    a, b = _dense(a), _dense(b)
    if a.shape != b.shape:
        return False

    missing_a, missing_b = np.isnan(a), np.isnan(b)
    if not np.array_equal(missing_a, missing_b):
        return False
    if missing_a.all():
        return True

    return bool(np.allclose(a[~missing_a], b[~missing_b], rtol=0, atol=_TOLERANCE))


def _equal_column(a: pd.Series, b: pd.Series) -> bool:
    """Compare two annotation columns, tolerating a change of storage type.

    Numeric-like columns are compared as floats, so the integer -> double cast
    that ``cast_nullable_columns()`` performs on the R side does not read as a
    value change; the cast itself is reported separately as a ``*_types``
    deviation. Everything else is compared as string, which absorbs
    categorical -> object without hiding a change of level.
    """
    if len(a) != len(b):
        return False

    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        left = pd.to_numeric(a, errors="coerce").to_numpy(dtype=float)
        right = pd.to_numeric(b, errors="coerce").to_numpy(dtype=float)
        return _equal_matrix(left, right)

    left = a.astype("object").where(a.notna(), None)
    right = b.astype("object").where(b.notna(), None)
    return [None if value is None else str(value) for value in left] == [
        None if value is None else str(value) for value in right
    ]


def _compare_frame(a: pd.DataFrame, b: pd.DataFrame, prefix: str, scope: str) -> list[Deviation]:
    """Compare two annotation frames column by column, positionally.

    Rows are compared by position rather than by label, matching how ``.X`` is
    compared. A change of index is reported once by the caller as a
    ``feature_names``/``obs_names`` deviation instead of once per column here.
    """
    deviations: list[Deviation] = []

    only_a, only_b = set(a.columns) - set(b.columns), set(b.columns) - set(a.columns)
    if only_a or only_b:
        deviations.append(
            Deviation(
                f"{prefix}_columns",
                scope,
                f"only in a: {{{','.join(sorted(only_a))}}}; only in b: {{{','.join(sorted(only_b))}}}",
            )
        )

    if len(a) != len(b):
        deviations.append(Deviation(f"{prefix}_values", scope, f"row count {len(a)} vs {len(b)}"))
        return deviations

    for column in [name for name in a.columns if name in set(b.columns)]:
        left, right = a[column].reset_index(drop=True), b[column].reset_index(drop=True)

        if str(left.dtype) != str(right.dtype):
            deviations.append(Deviation(f"{prefix}_types", scope, f"{column}: {left.dtype} vs {right.dtype}"))

        if not _equal_column(left, right):
            deviations.append(Deviation(f"{prefix}_values", scope, f"{column} differs"))

    return deviations


def _graph_edges(matrix: Any, names: Iterable[str]) -> dict[tuple[str, str], float]:
    """The ``.varp`` matrix as a name-keyed edge set.

    Keyed by feature name so that the graph can be
    compared independently of the order of the global feature axis.
    """
    names = list(names)
    coo = matrix.tocoo() if issparse(matrix) else None

    if coo is None:
        dense = _dense(matrix)
        rows, cols = np.nonzero(dense)
        values = dense[rows, cols]
    else:
        rows, cols, values = coo.row, coo.col, coo.data

    return {
        (names[int(row)], names[int(col)]): float(value)
        for row, col, value in zip(rows, cols, values, strict=True)
        if value != 0
    }


def _compare_modality(a: AnnData, b: AnnData, scope: str) -> list[Deviation]:
    deviations: list[Deviation] = []

    if a.shape != b.shape:
        # Everything below compares positionally, so a shape mismatch is
        # reported alone rather than amplified into a dozen deviations.
        return [Deviation("shape", scope, f"{a.shape} vs {b.shape}")]

    if set(a.layers) != set(b.layers):
        deviations.append(
            Deviation("layer_names", scope, f"a: [{','.join(sorted(a.layers))}]; b: [{','.join(sorted(b.layers))}]")
        )

    if list(a.var_names) != list(b.var_names):
        shared = sum(x == y for x, y in zip(a.var_names, b.var_names, strict=True))
        deviations.append(Deviation("feature_names", scope, f"{shared} of {a.n_vars} identical"))

    if list(a.obs_names) != list(b.obs_names):
        shared = sum(x == y for x, y in zip(a.obs_names, b.obs_names, strict=True))
        deviations.append(Deviation("obs_names", scope, f"{shared} of {a.n_obs} identical"))

    if not _equal_matrix(a.X, b.X):
        deviations.append(Deviation("values", scope, ".X differs"))

    deviations += _compare_frame(a.var, b.var, "feature_annotation", scope)
    deviations += _compare_frame(a.obs, b.obs, "obs_annotation", scope)

    return deviations


def compare_mudata(a: MuData, b: MuData) -> Report:
    """Compare two MuData objects against the round-trip conformance profile.

    Parameters
    ----------
    a
        The original object.
    b
        The object recovered from a round trip.

    Returns
    -------
    A :class:`Report` listing every property on which the two objects differ.
    """
    deviations: list[Deviation] = []

    names_a, names_b = list(a.mod), list(b.mod)
    if set(names_a) != set(names_b):
        deviations.append(
            Deviation(
                "sets",
                "",
                f"only in a: {{{','.join(sorted(set(names_a) - set(names_b)))}}}; "
                f"only in b: {{{','.join(sorted(set(names_b) - set(names_a)))}}}",
            )
        )
    elif names_a != names_b:
        deviations.append(Deviation("set_order", "", f"a: [{','.join(names_a)}]; b: [{','.join(names_b)}]"))

    if list(a.var_names) != list(b.var_names):
        deviations.append(Deviation("global_feature_names", "", f"a: {a.n_vars} names; b: {b.n_vars} names"))

    if list(a.obs_names) != list(b.obs_names):
        deviations.append(Deviation("global_obs_names", "", f"a: {a.n_obs} names; b: {b.n_obs} names"))

    deviations += _compare_frame(a.obs, b.obs, "global_obs_annotation", "")

    for name in [name for name in names_a if name in set(names_b)]:
        deviations += _compare_modality(a.mod[name], b.mod[name], name)

    keys_a, keys_b = set(a.varp), set(b.varp)
    if keys_a != keys_b:
        only_a = ",".join(sorted(keys_a - keys_b))
        only_b = ",".join(sorted(keys_b - keys_a))
        deviations.append(Deviation("feature_graph_keys", "", f"only in a: {{{only_a}}}; only in b: {{{only_b}}}"))

    for key in sorted(keys_a & keys_b):
        edges_a = _graph_edges(a.varp[key], a.var_names)
        edges_b = _graph_edges(b.varp[key], b.var_names)
        if edges_a != edges_b:
            shared = len(set(edges_a.items()) & set(edges_b.items()))
            deviations.append(
                Deviation("feature_graph", key, f"{len(edges_a)} vs {len(edges_b)} edges, {shared} identical")
            )

    return Report(deviations)
