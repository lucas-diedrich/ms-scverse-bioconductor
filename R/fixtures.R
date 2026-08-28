
## QFeatures fixtures for the round-trip integration tests.
##
## Every fixture is deterministic and reconstructible from its name alone, so a
## test can write it, round trip it through Python and rebuild the original for
## comparison without carrying an object across process boundaries.
##
## The datasets that ship with QFeatures cover the interesting cases without a
## download, which keeps the per-pull-request suite offline:
##
##   feat1   one assay, no links, globally unique row names, a global colData
##   feat2   three unlinked assays whose row names collide across assays
##   feat3   seven assays, a fan-in (psmsall <- psms1 + psms2), NA values in an
##           assay, and row names that collide across assays
##   ft_na   one assay carrying an NA pattern in the quantitative data
##
## `nullable_rowdata` is constructed rather than loaded, to exercise the
## NA-bearing integer/logical columns that `cast_nullable_columns()` handles.

library(QFeatures)
library(S4Vectors)

QFEATURES_DATASETS <- c("feat1", "feat2", "feat3", "feat4", "ft_na")

FIXTURES <- c(QFEATURES_DATASETS, "nullable_rowdata")


.fixture_dataset <- function(name) {
    env <- new.env(parent = emptyenv())
    utils::data(list = name, package = "QFeatures", envir = env)
    get(name, envir = env)
}


#' One assay whose rowData holds NA-bearing integer and logical columns.
#'
#' `MuData:::write_matrix()` routes those through the nullable-integer/boolean
#' encoding, whose mask anndata rejects, so `cast_nullable_columns()` casts them
#' to double. The fixture exists to pin that behaviour: the values and the
#' missingness must survive, the storage type must not.
.fixture_nullable_rowdata <- function() {
    qf <- .fixture_dataset("feat1")
    se <- qf[["psms"]]

    rd <- rowData(se)
    rd$Var[c(1L, 3L)] <- NA_integer_
    rd$flag <- rep(c(TRUE, FALSE, NA), length.out = nrow(rd))
    rowData(se) <- rd

    QFeatures(List(psms = se), colData = colData(qf))
}


#' Build a fixture by name.
#'
#' @param name One of `FIXTURES`.
#'
#' @return A `QFeatures` object.
fixture <- function(name) {
    if (identical(name, "nullable_rowdata"))
        return(.fixture_nullable_rowdata())

    if (!name %in% QFEATURES_DATASETS)
        stop("Unknown fixture '", name, "'. Available: '",
             paste(FIXTURES, collapse = "', '"), "'.")

    .fixture_dataset(name)
}
