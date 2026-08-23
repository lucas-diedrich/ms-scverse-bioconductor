
## Conformance comparison of two QFeatures objects.
##
## The R half of the round-trip conformance profile (ROADMAP.md section 4). The
## Python half lives in `ms_scverse_bioconductor.conformance` and reports the
## same record shape, so a test can assert on deviation kinds without caring
## which language produced them.
##
## A comparison returns a list of deviation records rather than a single
## boolean, because "lossless" is a profile and not a property: some deviations
## are known and accepted, and a test states which ones it tolerates. Reporting
## kinds separately is what keeps that list honest -- once a gap is closed the
## kind disappears and the test that still tolerates it fails.
##
##   compare_qfeatures(a, b)   -> list(equal = <logical>, deviations = <list>)

library(QFeatures)
library(MultiAssayExperiment)
library(S4Vectors)


## Deviation kinds. Shared with the Python comparator; only the subset that
## applies to a given language is ever emitted.
##
##   sets                        set of assay/modality names differs
##   set_order                   same names, different order
##   layer_names                 assayNames() of an assay differs
##   shape                       assay dimensions differ
##   feature_names               per-assay feature (row) names differ
##   obs_names                   per-assay observation (column) names differ
##   values                      quantitative values or missingness differ
##   feature_annotation_*        rowData columns / types / values
##   obs_annotation_*            per-assay colData columns / types / values
##   global_obs_annotation_*     QFeatures colData (MuData .obs)
##   sample_map                  MultiAssayExperiment sampleMap
##   link_parents                AssayLink@from, per assay
##   link_topology               feature-level parent -> child edges
##   link_fcol                   AssayLink@fcol


.deviation <- function(kind, scope, detail) {
    list(kind = kind, scope = scope, detail = detail)
}


#' Compare two numeric matrices, tolerating NA/NaN interchange.
#'
#' Missingness is compared as a pattern and values only where both are present.
#' `is.na()` is TRUE for NaN as well as NA, which is deliberate: MuData writes
#' an R `NA` as NaN, so the distinction cannot survive an .h5mu round trip and
#' is listed as a known loss.
.equal_matrix <- function(x, y, tolerance = 1e-8) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (!identical(dim(x), dim(y))) return(FALSE)

    nx <- is.na(x)
    ny <- is.na(y)
    if (!identical(unname(nx), unname(ny))) return(FALSE)
    if (all(nx)) return(TRUE)

    isTRUE(all.equal(unname(x[!nx]), unname(y[!ny]), tolerance = tolerance))
}


#' Compare two annotation columns, tolerating a change of storage type.
#'
#' Numeric-like columns are compared as doubles, so the integer -> double cast
#' that `cast_nullable_columns()` performs does not read as a value change; the
#' cast itself is reported separately as a `*_types` deviation. Everything else
#' is compared as character, which absorbs factor -> character without hiding a
#' change of level.
.equal_column <- function(a, b, tolerance = 1e-8) {
    numeric_like <- function(v) is.numeric(v) || is.logical(v)

    if (numeric_like(a) && numeric_like(b)) {
        a <- as.double(a)
        b <- as.double(b)
        na <- is.na(a)
        nb <- is.na(b)
        if (!identical(na, nb)) return(FALSE)
        if (all(na)) return(TRUE)
        return(isTRUE(all.equal(a[!na], b[!nb], tolerance = tolerance)))
    }

    a <- as.character(a)
    b <- as.character(b)
    identical(is.na(a), is.na(b)) && identical(a[!is.na(a)], b[!is.na(b)])
}


#' Compare two DataFrames column by column, positionally.
#'
#' Rows are compared by position rather than by name, matching how the assay
#' matrices are compared. A change of row names is reported once, by the caller,
#' as a `feature_names`/`obs_names` deviation instead of once per column here.
.compare_dataframe <- function(x, y, prefix, scope) {
    dev <- list()

    cx <- colnames(x)
    cy <- colnames(y)
    if (!setequal(cx, cy))
        dev <- c(dev, list(.deviation(
            paste0(prefix, "_columns"), scope,
            sprintf("only in a: {%s}; only in b: {%s}",
                    paste(setdiff(cx, cy), collapse = ","),
                    paste(setdiff(cy, cx), collapse = ",")))))

    if (nrow(x) != nrow(y)) {
        dev <- c(dev, list(.deviation(
            paste0(prefix, "_values"), scope,
            sprintf("row count %d vs %d", nrow(x), nrow(y)))))
        return(dev)
    }

    for (column in intersect(cx, cy)) {
        a <- x[[column]]
        b <- y[[column]]

        class_a <- class(a)[1]
        class_b <- class(b)[1]
        if (!identical(class_a, class_b))
            dev <- c(dev, list(.deviation(
                paste0(prefix, "_types"), scope,
                sprintf("%s: %s vs %s", column, class_a, class_b))))

        if (!.equal_column(a, b))
            dev <- c(dev, list(.deviation(
                paste0(prefix, "_values"), scope,
                sprintf("%s differs", column))))
    }

    dev
}


#' Canonical parent -> child edge set of one assay's `AssayLink`.
#'
#' `mcols(hits)$names_from`/`names_to` are the authoritative endpoints; the
#' integer indices and `nLnode`/`nRnode` are stale on any subset object, since
#' `QFeatures:::.pruneHits()` drops hits by row name without renumbering them.
#' Comparing the sorted set of `parent|from->to` strings therefore also makes
#' the comparison insensitive to the order hits happen to be stored in.
.link_edges <- function(object, name) {
    al <- assayLink(object, name)
    hits <- if (is(al@hits, "List")) as.list(al@hits) else list(al@hits)

    edges <- character(0)
    for (k in seq_along(hits)) {
        parent <- al@from[k]
        h <- hits[[k]]
        if (length(parent) == 0L || is.na(parent) || length(h) == 0L) next

        md <- S4Vectors::mcols(h)
        edges <- c(edges, sprintf("%s|%s->%s", parent,
                                  md$names_from, md$names_to))
    }

    sort(unique(edges))
}


.sample_map_rows <- function(object) {
    sm <- sampleMap(object)
    sort(sprintf("%s|%s|%s",
                 as.character(sm$assay),
                 as.character(sm$primary),
                 as.character(sm$colname)))
}


#' Compare two `QFeatures` objects against the round-trip conformance profile.
#'
#' @param a The original object.
#' @param b The object recovered from a round trip.
#'
#' @return A list with `equal` and a list of deviation records, each holding a
#'     `kind`, the `scope` it applies to (an assay name, or `""` for the object
#'     as a whole) and a human readable `detail`.
compare_qfeatures <- function(a, b) {
    dev <- list()

    names_a <- names(a)
    names_b <- names(b)

    if (!setequal(names_a, names_b))
        dev <- c(dev, list(.deviation(
            "sets", "",
            sprintf("only in a: {%s}; only in b: {%s}",
                    paste(setdiff(names_a, names_b), collapse = ","),
                    paste(setdiff(names_b, names_a), collapse = ",")))))
    else if (!identical(names_a, names_b))
        dev <- c(dev, list(.deviation(
            "set_order", "",
            sprintf("a: [%s]; b: [%s]",
                    paste(names_a, collapse = ","),
                    paste(names_b, collapse = ",")))))

    dev <- c(dev, .compare_dataframe(colData(a), colData(b),
                                     "global_obs_annotation", ""))

    map_a <- .sample_map_rows(a)
    map_b <- .sample_map_rows(b)
    if (!identical(map_a, map_b))
        dev <- c(dev, list(.deviation(
            "sample_map", "",
            sprintf("%d vs %d rows, %d shared",
                    length(map_a), length(map_b),
                    length(intersect(map_a, map_b))))))

    for (name in intersect(names_a, names_b)) {
        se_a <- a[[name]]
        se_b <- b[[name]]

        if (!identical(dim(se_a), dim(se_b))) {
            ## Everything below compares positionally, so a shape mismatch is
            ## reported alone rather than amplified into a dozen deviations.
            dev <- c(dev, list(.deviation(
                "shape", name,
                sprintf("%s vs %s",
                        paste(dim(se_a), collapse = "x"),
                        paste(dim(se_b), collapse = "x")))))
            next
        }

        ## Quoted, because the interesting case is a name that became the empty
        ## string: MuData stores the primary matrix as .X, which carries no name.
        if (!setequal(assayNames(se_a), assayNames(se_b)))
            dev <- c(dev, list(.deviation(
                "layer_names", name,
                sprintf("a: [%s]; b: [%s]",
                        paste(sprintf("'%s'", assayNames(se_a)), collapse = ","),
                        paste(sprintf("'%s'", assayNames(se_b)), collapse = ",")))))

        if (!identical(rownames(se_a), rownames(se_b)))
            dev <- c(dev, list(.deviation(
                "feature_names", name,
                sprintf("%d of %d identical",
                        sum(rownames(se_a) == rownames(se_b)),
                        nrow(se_a)))))

        if (!identical(colnames(se_a), colnames(se_b)))
            dev <- c(dev, list(.deviation(
                "obs_names", name,
                sprintf("%d of %d identical",
                        sum(colnames(se_a) == colnames(se_b)),
                        ncol(se_a)))))

        if (!.equal_matrix(assay(se_a), assay(se_b)))
            dev <- c(dev, list(.deviation("values", name, "assay differs")))

        dev <- c(dev, .compare_dataframe(rowData(se_a), rowData(se_b),
                                         "feature_annotation", name))
        dev <- c(dev, .compare_dataframe(colData(se_a), colData(se_b),
                                         "obs_annotation", name))

        parents_a <- sort(setdiff(assayLink(a, name)@from, NA_character_))
        parents_b <- sort(setdiff(assayLink(b, name)@from, NA_character_))
        if (!identical(parents_a, parents_b))
            dev <- c(dev, list(.deviation(
                "link_parents", name,
                sprintf("a: [%s]; b: [%s]",
                        paste(parents_a, collapse = ","),
                        paste(parents_b, collapse = ",")))))

        edges_a <- .link_edges(a, name)
        edges_b <- .link_edges(b, name)
        if (!identical(edges_a, edges_b))
            dev <- c(dev, list(.deviation(
                "link_topology", name,
                sprintf("%d vs %d edges, %d shared",
                        length(edges_a), length(edges_b),
                        length(intersect(edges_a, edges_b))))))

        fcol_a <- sort(setdiff(assayLink(a, name)@fcol, NA_character_))
        fcol_b <- sort(setdiff(assayLink(b, name)@fcol, NA_character_))
        if (!identical(fcol_a, fcol_b))
            dev <- c(dev, list(.deviation(
                "link_fcol", name,
                sprintf("a: [%s]; b: [%s]",
                        paste(fcol_a, collapse = ","),
                        paste(fcol_b, collapse = ",")))))
    }

    list(equal = length(dev) == 0L, deviations = dev)
}
