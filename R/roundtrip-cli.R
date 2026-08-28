
## Command line entry point for the R leg of the round-trip integration tests.
##
## The R and Python stacks talk to each other only through .h5mu files, never
## in-process. That is the path users actually take, it keeps rhdf5 and h5py in
## separate processes, and it means a failure is attributable to one leg. This
## script is the whole R-side surface the pytest harness drives, and it is
## equally usable by hand or through `docker run`.
##
##   Rscript R/roundtrip-cli.R check
##   Rscript R/roundtrip-cli.R list-fixtures
##   Rscript R/roundtrip-cli.R fixture    --name feat3   --out r0.h5mu
##   Rscript R/roundtrip-cli.R read-write --in a.h5mu    --out b.h5mu
##   Rscript R/roundtrip-cli.R compare    --fixture feat3 --in r1.h5mu
##   Rscript R/roundtrip-cli.R selftest
##
## `compare` writes its report to stdout as JSON, or to `--report <path>`. It
## exits 0 whenever the comparison ran, whatever the verdict: the report carries
## the verdict, and the exit code distinguishes "compared" from "R failed".

.script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg))
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    getwd()
}


#' Parse `subcommand --key value --flag` without a dependency.
#'
#' optparse is not part of a Bioconductor base image and the grammar here is one
#' subcommand plus long options, so hand-parsing avoids an install for nothing.
.parse_args <- function(argv) {
    if (!length(argv))
        stop("Missing subcommand. See the header of this file for usage.")

    parsed <- list(command = argv[1])
    rest <- argv[-1]

    i <- 1L
    while (i <= length(rest)) {
        token <- rest[i]
        if (!startsWith(token, "--"))
            stop("Expected an option starting with '--', got '", token, "'.")

        key <- sub("^--", "", token)
        has_value <- i + 1L <= length(rest) && !startsWith(rest[i + 1L], "--")
        if (has_value) {
            parsed[[key]] <- rest[i + 1L]
            i <- i + 2L
        } else {
            parsed[[key]] <- TRUE
            i <- i + 1L
        }
    }

    parsed
}


.require_option <- function(args, key) {
    value <- args[[key]]
    if (is.null(value) || isTRUE(value))
        stop("Option '--", key, "' is required and needs a value.")
    value
}


.emit <- function(payload, path = NULL) {
    json <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null",
                             digits = NA, pretty = TRUE)
    if (is.null(path)) cat(json, "\n", sep = "") else writeLines(json, path)
    invisible(NULL)
}


## `check` must run before anything is sourced, so that a missing dependency is
## reported as a readable payload instead of an error from a library() call.
.cmd_check <- function(args) {
    required <- c("QFeatures", "MuData", "MultiAssayExperiment", "rhdf5",
                  "Matrix", "S4Vectors", "jsonlite")

    versions <- vapply(required, function(p) {
        tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
    }, character(1))

    missing <- names(versions)[is.na(versions)]
    payload <- list(ok = length(missing) == 0L,
                    r_version = paste(R.version$major, R.version$minor, sep = "."),
                    packages = as.list(versions),
                    missing = as.list(missing))

    .emit(payload, if (is.character(args$report)) args$report else NULL)
    if (!payload$ok) quit(status = 1L)
}


.cmd_list_fixtures <- function(args) {
    .emit(list(fixtures = as.list(FIXTURES)),
          if (is.character(args$report)) args$report else NULL)
}


.cmd_fixture <- function(args) {
    name <- .require_option(args, "name")
    out <- .require_option(args, "out")

    object <- fixture(name)
    writeLinkH5MU(object, out, overwrite = TRUE)
    message("Wrote ", out, " from fixture '", name, "'.")
}


.cmd_read_write <- function(args) {
    ## `in` is a reserved word, so it is only reachable through [[.
    input <- .require_option(args, "in")
    out <- .require_option(args, "out")

    object <- readLinkH5MU(input)
    writeLinkH5MU(object, out, overwrite = TRUE)
    message("Read ", input, " into a QFeatures object with ",
            length(object), " assays and wrote ", out, ".")
}


.cmd_selftest <- function(args) {
    report <- compare_selftest()
    .emit(report, if (is.character(args$report)) args$report else NULL)
    if (!report$ok) quit(status = 1L)
}


.cmd_compare <- function(args) {
    name <- .require_option(args, "fixture")
    input <- .require_option(args, "in")

    report <- compare_qfeatures(fixture(name), readLinkH5MU(input))
    .emit(report, if (is.character(args$report)) args$report else NULL)
}


main <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- .parse_args(argv)

    if (identical(args$command, "check")) {
        if (!requireNamespace("jsonlite", quietly = TRUE))
            stop("jsonlite is required to report the result of 'check'.")
        return(.cmd_check(args))
    }

    here <- .script_dir()
    suppressPackageStartupMessages({
        source(file.path(here, "io.R"))
        source(file.path(here, "compare.R"))
        source(file.path(here, "fixtures.R"))
        source(file.path(here, "compare-selftest.R"))
    })

    switch(args$command,
           "list-fixtures" = .cmd_list_fixtures(args),
           "fixture"       = .cmd_fixture(args),
           "read-write"    = .cmd_read_write(args),
           "compare"       = .cmd_compare(args),
           "selftest"      = .cmd_selftest(args),
           stop("Unknown subcommand '", args$command, "'."))
}


if (!interactive()) main()
