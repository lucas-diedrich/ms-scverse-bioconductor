
## Self-test for compare_qfeatures().
##
## The round-trip profiles in tests/integration assert that a comparison reports
## exactly a given set of deviation kinds. A comparator that reported nothing
## would make every one of them pass vacuously, so each kind it can emit is
## provoked here by a known mutation and checked.
##
## This is the R counterpart of tests/test_conformance.py and is driven through
## `Rscript R/roundtrip-cli.R selftest`, because testthat is not part of a
## Bioconductor base image and one JSON payload is all the harness needs.

library(QFeatures)
library(MultiAssayExperiment)
library(S4Vectors)


#' Replace one assay without triggering the QFeatures validity check.
#'
#' The mutations below deliberately produce objects whose `AssayLinks` no longer
#' match their assays, which is exactly what a broken round trip looks like, so
#' the slot is assigned directly rather than through `[[<-`.
.replace_assay <- function(object, name, se) {
    assays <- as.list(experiments(object))
    assays[[name]] <- se
    object@ExperimentList <- ExperimentList(assays)
    object
}


.mutate_value <- function(object) {
    se <- object[["psms"]]
    assay(se)[1L, 1L] <- assay(se)[1L, 1L] + 1
    .replace_assay(object, "psms", se)
}


.mutate_missing <- function(object) {
    se <- object[["psms"]]
    assay(se)[2L, 1L] <- NA_real_
    .replace_assay(object, "psms", se)
}


.mutate_feature_name <- function(object) {
    se <- object[["psms"]]
    rownames(se)[1L] <- "renamed"
    .replace_assay(object, "psms", se)
}


.mutate_annotation_column <- function(object) {
    se <- object[["psms"]]
    rowData(se)$extra <- seq_len(nrow(se))
    .replace_assay(object, "psms", se)
}


.mutate_annotation_type <- function(object) {
    se <- object[["psms"]]
    ## Same values, different storage type: a type deviation and nothing else.
    rowData(se)$Var <- as.double(rowData(se)$Var)
    .replace_assay(object, "psms", se)
}


.mutate_observation_name <- function(object) {
    se <- object[["psms"]]
    colnames(se)[1L] <- "renamed"
    .replace_assay(object, "psms", se)
}


.mutate_shape <- function(object) {
    .replace_assay(object, "psms", object[["psms"]][seq_len(3L), ])
}


.mutate_link_edge <- function(object) {
    al <- object@assayLinks[["peptides"]]
    al@hits <- al@hits[-1L]
    object@assayLinks[["peptides"]] <- al
    object
}


.mutate_link_fcol <- function(object) {
    al <- object@assayLinks[["peptides"]]
    al@fcol <- "changed"
    object@assayLinks[["peptides"]] <- al
    object
}


.mutate_link_parent <- function(object) {
    al <- object@assayLinks[["peptides"]]
    al@from <- "psms1"
    object@assayLinks[["peptides"]] <- al
    object
}


#' Cases the self-test runs. Each names the fixture to start from and either a
#' mutation to apply or a second fixture to compare against.
.selftest_cases <- function() {
    list(
        list(name = "identical",         a = "feat1", expected = character(0)),
        list(name = "changed_value",     a = "feat1", expected = "values",
             mutate = .mutate_value),
        list(name = "introduced_na",     a = "feat1", expected = "values",
             mutate = .mutate_missing),
        list(name = "renamed_feature",   a = "feat1", expected = "feature_names",
             mutate = .mutate_feature_name),
        ## Renaming a column of an assay leaves the sampleMap slot untouched,
        ## so this provokes obs_names alone; sample_map is provoked below.
        list(name = "renamed_obs",       a = "feat1", expected = "obs_names",
             mutate = .mutate_observation_name),
        list(name = "added_column",      a = "feat1",
             expected = "feature_annotation_columns",
             mutate = .mutate_annotation_column),
        list(name = "changed_type",      a = "feat1",
             expected = "feature_annotation_types",
             mutate = .mutate_annotation_type),
        list(name = "changed_shape",     a = "feat1", expected = "shape",
             mutate = .mutate_shape),
        list(name = "dropped_link_edge", a = "feat3", expected = "link_topology",
             mutate = .mutate_link_edge),
        list(name = "changed_fcol",      a = "feat3", expected = "link_fcol",
             mutate = .mutate_link_fcol),
        list(name = "changed_parent",    a = "feat3",
             expected = c("link_parents", "link_topology"),
             mutate = .mutate_link_parent),
        ## Two unrelated objects: everything that is compared at the top level
        ## differs, which is what pins the global kinds.
        list(name = "different_object",  a = "feat1", b = "feat3",
             expected = c("sets", "sample_map",
                          "global_obs_annotation_columns",
                          "global_obs_annotation_values"))
    )
}


#' Run every self-test case.
#'
#' @return A list with `ok` and one record per case, holding the expected and
#'     the observed deviation kinds.
compare_selftest <- function() {
    results <- lapply(.selftest_cases(), function(case) {
        a <- fixture(case$a)
        b <- if (!is.null(case$b)) fixture(case$b)
             else if (!is.null(case$mutate)) case$mutate(fixture(case$a))
             else fixture(case$a)

        report <- compare_qfeatures(a, b)
        observed <- sort(unique(vapply(report$deviations, `[[`, character(1), "kind")))
        expected <- sort(case$expected)

        list(name = case$name,
             ok = identical(observed, expected),
             expected = as.list(expected),
             observed = as.list(observed))
    })

    list(ok = all(vapply(results, `[[`, logical(1), "ok")), cases = results)
}
