
## Convert between .h5mu files and QFeatures objects, in both directions.
##
## MuData::readH5MU()/writeH5MU() cover the MultiAssayExperiment skeleton
## (assays, colData, sampleMap, uns). Neither touches the *global* .varp group,
## because MultiAssayExperiment has nowhere to put a pairwise feature matrix.
## That group is where mulink keeps the feature graph, so it is read and written
## separately, against a QFeatures object, whose AssayLinks can represent it.
##
##   qfeatures_read_h5mu()   .h5mu -> QFeatures
##   writeQFeaturesH5MU()    QFeatures -> .h5mu

library(QFeatures)
library(MuData)
library(MultiAssayExperiment)
library(rhdf5)
library(Matrix)
library(S4Vectors)

## Internals of MuData that a patched readH5MU() would call directly.
.h5autoclose <- MuData:::h5autoclose
.read_matrix <- MuData:::read_matrix
.read_with_index <- MuData:::read_with_index
.write_matrix <- MuData:::write_matrix
.write_data_frame <- MuData:::write_data_frame


#' Global feature names of an .h5mu file, in file order.
.var_names <- function(h5) {
    rownames(.read_with_index(.h5autoclose(h5 & "var")))
}


#' Read the global `.varp` group of an .h5mu file.
#'
#' Mirrors the `varp` branch of `MuData:::read_modality()`, which reads a
#' modality's pairwise feature matrices into `rowPair()`. Here the matrices are
#' returned instead, since they span all modalities and belong to no single one.
#'
#' `read_matrix()` transposes to maintain directionality as observation/feature directionality is inversed between
#' scverse and Bioconductor functions
#'
#' @param file Path to the .h5mu file.
#' @param keys Names to read from `.varp`. All of them if `NULL`.
#'
#' @return A named list of sparse matrices, with the global feature names as
#'     dimnames.
read_varp <- function(file, keys = NULL) {
    h5 <- H5Fopen(file, flags = "H5F_ACC_RDONLY", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    if (!H5Lexists(h5, "varp")) return(list())

    var_names <- .var_names(h5)
    available <- h5ls(.h5autoclose(h5 & "varp"), recursive = FALSE)$name
    if (is.null(keys)) keys <- available

    missing_keys <- setdiff(keys, available)
    if (length(missing_keys))
        warning("No '", paste(missing_keys, collapse = "', '"), "' in .varp.")

    matrices <- lapply(setNames(nm = intersect(keys, available)), function(key) {
        m <- .read_matrix(.h5autoclose(h5 & paste("varp", key, sep = "/")))
        if (!is(m, "dsparseMatrix")) {
            warning("Pairwise varp matrix ", key, " is not a sparse matrix. ",
                    "Only sparse matrices are currently supported, skipping...")
            return(NULL)
        }
        m <- t(m)
        dimnames(m) <- list(var_names, var_names)
        m
    })

    matrices[!vapply(matrices, is.null, logical(1))]
}


#' Build a `Hits` object from an adjacency submatrix.
#'
#' Rows are parent features, columns are child features. `names_from`/`names_to`
#' are required by QFeatures' validity checks (`.checkLinksInHits`), which verify
#' that every linked feature exists in the corresponding assay.
#'
#' `sort.by.query = TRUE` yields a `SortedByQueryHits`, which is what
#' `findMatches()` and therefore the rest of QFeatures produces.
hits_from_adjacency <- function(adj) {
    nz <- Matrix::summary(as(adj, "TsparseMatrix"))
    nz <- nz[order(nz$i, nz$j), , drop = FALSE]

    hits <- S4Vectors::Hits(from = nz$i,
                            to = nz$j,
                            nLnode = nrow(adj),
                            nRnode = ncol(adj),
                            sort.by.query = TRUE)

    S4Vectors::mcols(hits)$names_from <- rownames(adj)[nz$i]
    S4Vectors::mcols(hits)$names_to <- colnames(adj)[nz$j]
    hits
}


#' Turn a feature-mapping matrix into one `AssayLink` per assay.
#'
#' An edge u -> v in `.varp` means v's assay was derived from u's, so edges are
#' grouped by their *target* modality: that is the child, and the modalities the
#' edges come from are its parents. Grouping matters — `AssayLinks` is keyed by
#' child name, so adding links pair-by-pair would keep only the last parent of a
#' fan-in.
assay_links_from_feature_mapping <- function(experiments,
                                             feature_mapping,
                                             fcol = NA_character_) {
    mod_names <- names(experiments)
    var_names <- rownames(feature_mapping)

    ## The mulink convention requires globally unique feature names, so matching
    ## each modality's rownames against the global index is unambiguous.
    positions <- lapply(experiments, function(se) match(rownames(se), var_names))
    unmatched <- vapply(positions, anyNA, logical(1))
    if (any(unmatched))
        stop("Features of modality/modalities '",
             paste(mod_names[unmatched], collapse = "', '"),
             "' are absent from the global .var index.")

    has_edges <- function(adj) length(adj@x) > 0L && any(adj@x != 0)

    edges <- list()
    for (from in mod_names) {
        for (to in setdiff(mod_names, from)) {
            adj <- feature_mapping[positions[[from]], positions[[to]], drop = FALSE]
            if (!has_edges(adj)) next

            dimnames(adj) <- list(rownames(experiments[[from]]),
                                  rownames(experiments[[to]]))
            edges[[to]] <- c(edges[[to]], setNames(list(adj), from))
        }
    }

    within_modality <- vapply(mod_names, function(m) {
        has_edges(feature_mapping[positions[[m]], positions[[m]], drop = FALSE])
    }, logical(1))
    if (any(within_modality))
        warning("Ignoring edges within modality/modalities '",
                paste(mod_names[within_modality], collapse = "', '"),
                "'; an assay cannot be linked to itself.")

    ## One AssayLink per assay, in assay order: QFeatures' validity checks
    ## compare names(assayLinks) to names(object) with identical().
    AssayLinks(lapply(mod_names, function(to) {
        parents <- edges[[to]]
        if (is.null(parents)) return(AssayLink(name = to))

        hits <- lapply(parents, hits_from_adjacency)
        if (length(hits) > 1L) {
            hits <- S4Vectors::List(hits)
            names(hits) <- names(parents)
        } else {
            hits <- hits[[1L]]
        }

        AssayLink(name = to,
                  from = names(parents),
                  fcol = rep(fcol, length(parents)),
                  hits = hits)
    }))
}


#' Read an .h5mu file and create a `QFeatures` object.
#'
#' @param path Path to the .h5mu file.
#' @param feature_mapping_key Key of the feature graph in the global `.varp`.
#' @param backed Passed to `MuData::readH5MU()`.
#'
#' @return A `QFeatures` object. Assays and their links come from the file; if
#'     the `.varp` key is absent the assays are returned unlinked.
qfeatures_read_h5mu <- function(path,
                                feature_mapping_key = "feature_mapping",
                                backed = FALSE) {
    mae <- MuData::readH5MU(path, backed = backed)
    feature_mapping <- read_varp(path, keys = feature_mapping_key)[[feature_mapping_key]]

    assay_links <- if (is.null(feature_mapping)) {
        AssayLinks(names = names(experiments(mae)))
    } else {
        ## fcol normally names the rowData column a link was derived from. The
        ## link comes from .varp here, so record that key instead; the rowData
        ## column that produced it upstream lives in .uns (see ROADMAP.md).
        assay_links_from_feature_mapping(experiments(mae),
                                         feature_mapping,
                                         fcol = feature_mapping_key)
    }

    ## Promotion mirrors the QFeatures() constructor, which builds a
    ## MultiAssayExperiment and copies its slots across.
    qf <- new("QFeatures",
              ExperimentList = experiments(mae),
              colData = colData(mae),
              sampleMap = sampleMap(mae),
              metadata = metadata(mae),
              assayLinks = assay_links)

    validObject(qf)
    qf
}


#' Global feature axis of a `QFeatures` object.
#'
#' `MuData::writeH5MU()` builds the global `.var` index by `rbind`-ing one
#' zero-column data.frame per modality in assay order, so the axis is the assays'
#' row names concatenated in that order. This is the axis a `.varp` matrix has to
#' be aligned to.
#'
#' QFeatures enforces row-name uniqueness only *within* an assay
#' (`.unique_row_names`), while mudata keys the global axis by name and mulink
#' queries resolve it with `var_names.get_indexer()`, which fails on duplicates.
#' Colliding names are therefore prefixed with the assay name and the original ID
#' is kept in `.var` (see ROADMAP.md). Prefixing is all-or-nothing rather than
#' per-feature, so that the axis stays predictable.
#'
#' @param object A `QFeatures` object.
#' @param prefix When to prefix feature names with the assay name: `"collision"`
#'     only if the plain names are not globally unique, `"always"`, or `"never"`.
#' @param sep Separator between the assay name and the feature name.
#'
#' @return A list of the global feature `key`s and, parallel to them, the `assay`
#'     and original `id` each came from.
feature_index <- function(object,
                          prefix = c("collision", "always", "never"),
                          sep = ":") {
    prefix <- match.arg(prefix)

    ids <- lapply(names(object), function(nm) rownames(object[[nm]]))
    names(ids) <- names(object)
    assay <- rep(names(ids), lengths(ids))
    ids <- unlist(ids, use.names = FALSE)

    key <- switch(prefix,
                  always = paste0(assay, sep, ids),
                  never = ids,
                  collision = if (anyDuplicated(ids)) paste0(assay, sep, ids) else ids)

    if (anyDuplicated(key))
        stop("Global feature names are not unique, even after prefixing with ",
             "the assay name. Duplicates: '",
             paste(unique(key[duplicated(key)]), collapse = "', '"), "'.")

    list(key = key, assay = assay, id = ids)
}


#' Build the feature-mapping matrix from a `QFeatures` object's `AssayLinks`.
#'
#' Inverse of `assay_links_from_feature_mapping()`. Every `Hits` object in
#' `AssayLinks` is already a sparse edge list, so the only work is relabelling its
#' endpoints from per-assay row names to positions on the global feature axis.
#'
#' Edges are directed parent -> child: `A[u, v] != 0` means the feature `v` was
#' derived from the feature `u`. Only direct edges are stored, never the
#' transitive closure — mulink traverses multiple hops at query time, and a
#' closed graph would not survive the round trip, since transitive reduction is
#' not unique (see rfc.md and ROADMAP.md).
#'
#' @param object A `QFeatures` object.
#' @param index The global feature axis, from `feature_index()`.
#'
#' @return A `p x p` sparse matrix over the global feature axis.
feature_mapping_from_assay_links <- function(object,
                                             index = feature_index(object)) {
    p <- length(index$key)

    ## Assay name -> global positions, named by that assay's row names, so that a
    ## feature name resolves to a row/column of the matrix by lookup.
    positions <- split(seq_len(p), factor(index$assay, levels = names(object)))
    positions <- mapply(function(pos, nm) setNames(pos, rownames(object[[nm]])),
                        positions, names(positions), SIMPLIFY = FALSE)

    edges <- lapply(names(object), function(to) {
        al <- assayLink(object, to)

        ## @hits is either a single Hits or a List of Hits parallel to @from.
        ## Both shapes occur even with a single parent: the adjacency-matrix
        ## branch of QFeatures:::.create_assay_link() wraps its result with
        ## `if (length(hits) > 1)` while `hits` is still a Hits object, so the
        ## test counts hits rather than parents and yields a length-1 HitsList.
        hits <- if (is(al@hits, "List")) as.list(al@hits) else list(al@hits)

        do.call(rbind, lapply(seq_along(hits), function(k) {
            from <- al@from[k]
            h <- hits[[k]]
            if (is.na(from) || length(h) == 0L) return(NULL)  # root assay
            md <- mcols(h)

            ## names_from/names_to are the authoritative endpoints.
            ## QFeatures:::.pruneHits() drops hits by row name when the object is
            ## subset but never renumbers the integer from/to indices, nor
            ## nLnode/nRnode, so those are stale on any subset object.
            data.frame(i = positions[[from]][md$names_from],
                       j = positions[[to]][md$names_to],
                       x = if ("weigths" %in% colnames(md)) as.numeric(md$weigths)
                           else rep(1, length(h)))
        }))
    })

    edges <- do.call(rbind, edges)
    if (is.null(edges))
        edges <- data.frame(i = integer(0), j = integer(0), x = numeric(0))

    sparseMatrix(i = edges$i, j = edges$j, x = edges$x,
                 dims = c(p, p), dimnames = list(index$key, index$key),
                 repr = "R")
}


#' Rewrite the global `.var` index of an .h5mu file.
#'
#' `MuData::writeH5MU()` builds the global index as
#' `do.call(rbind, vars)` over one *zero-column* data.frame per modality
#' (`write_h5mu.R:184`). `rbind` on zero-column data.frames returns zero rows, so
#' the group it writes is always empty, whatever the object contained. Verified
#' against MuData 1.14.0.
#'
#' The index is not optional here: mudata keys the global axis by name, and both
#' `read_varp()` and mulink resolve `.varp` against it. So it is written again,
#' correctly, after `writeH5MU()` has returned.
#'
#' @param file Path to the .h5mu file.
#' @param keys The global feature names, in file order.
write_var_index <- function(file, keys) {
    h5 <- H5Fopen(file, flags = "H5F_ACC_RDWR", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    if (H5Lexists(h5, "var")) H5Ldelete(h5, "var")
    .write_data_frame(h5, "var", data.frame(row.names = keys))
    invisible(NULL)
}


#' Append a matrix to the global `.varp` group of an .h5mu file.
#'
#' Counterpart of `read_varp()`. `MuData::writeH5MU()` writes no global `.varp`,
#' so the group is added afterwards. By then the file is closed and finalized
#' (`finalize_mudata()`); reopening it read/write preserves the 512-byte user
#' block that holds the MuData signature.
#'
#' `write_matrix()` stores a `dgRMatrix`'s row-compressed arrays under the
#' `csr_matrix` encoding, which is exactly how Python reads them back, so
#' `mdata.varp[key][u, v]` in Python is `mat[u, v]` here. `read_varp()` inverts
#' this: `read_sparse_matrix()` reads a `csr_matrix` as its transpose, and
#' `read_varp()` transposes again.
#'
#' @param file Path to the .h5mu file.
#' @param key Name to store the matrix under in `.varp`.
#' @param mat A sparse matrix over the global feature axis.
write_varp <- function(file, key, mat) {
    h5 <- H5Fopen(file, flags = "H5F_ACC_RDWR", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    if (H5Lexists(h5, "varp")) {
        grp <- H5Gopen(h5, "varp")
        if (H5Lexists(grp, key)) H5Ldelete(grp, key)
    } else {
        grp <- H5Gcreate(h5, "varp")
    }
    ## after = FALSE so the group is closed before the file.
    on.exit(H5Gclose(grp), add = TRUE, after = FALSE)

    .write_matrix(grp, key, as(mat, "RsparseMatrix"), needTranspose = TRUE)
    invisible(NULL)
}


#' Write a `QFeatures` object to an .h5mu file.
#'
#' Inverse of `qfeatures_read_h5mu()`. `MuData::writeH5MU()` covers the
#' MultiAssayExperiment skeleton; the feature graph held in `AssayLinks` has no
#' place there, so it is written separately to the global `.varp` as a `p x p`
#' sparse adjacency matrix over the global `.var` index.
#'
#' Two parts of the specification are deliberately not implemented here and are
#' tracked separately in ROADMAP.md:
#'
#'   * The modality/layer split. Every assay is written as its own modality, so a
#'     one-to-one transformation such as `logTransform()` becomes a second
#'     modality plus an edge, rather than a layer of its parent.
#'   * `uns["mulink"]["assays"]`, which records `AssayLink@fcol` and which layer
#'     an inter-modality edge starts from. `fcol` therefore does not survive this
#'     write, and `qfeatures_read_h5mu()` substitutes the `.varp` key for it.
#'
#' Matrix-valued `rowData` columns are dropped with a warning; see below.
#'
#' @param object A `QFeatures` object.
#' @param path Path of the .h5mu file to create.
#' @param feature_mapping_key Key to store the feature graph under in `.varp`.
#' @param prefix,sep Passed to `feature_index()`.
#' @param overwrite Whether to replace `path` if it already exists.
#'
#' @return `path`, invisibly.
writeQFeaturesH5MU <- function(object, path,
                               feature_mapping_key = "feature_mapping",
                               prefix = c("collision", "always", "never"),
                               sep = ":",
                               overwrite = FALSE) {
    if (!is(object, "QFeatures"))
        stop("'object' must be a QFeatures object.")
    if (file.exists(path)) {
        ## MuData::writeH5MU() goes through H5Fcreate() and cannot truncate.
        if (!overwrite)
            stop("'", path, "' already exists; pass overwrite = TRUE to replace it.")
        unlink(path)
    }

    index <- feature_index(object, prefix = prefix, sep = sep)
    ## Built before any renaming, since the AssayLinks are keyed on the row names
    ## the object currently carries.
    feature_mapping <- feature_mapping_from_assay_links(object, index)

    keys <- split(index$key, factor(index$assay, levels = names(object)))
    dropped <- character(0)

    prepared <- mapply(function(se, name, key) {
        ## MuData cannot encode a matrix-valued rowData column: writeH5AD()
        ## stores it as a nested sparse-matrix group inside `var`, which
        ## read_dataframe() then cannot coerce back to a DataFrame column, so the
        ## file it produces is unreadable. Such a column is an adjacency matrix
        ## (see QFeatures:::.create_assay_link()) -- the same graph being written
        ## to .varp -- so it is dropped rather than left to break the file.
        rd <- rowData(se)
        is_matrix <- vapply(rd, function(x) is.matrix(x) || is(x, "Matrix"),
                            logical(1))
        if (any(is_matrix)) {
            dropped <<- c(dropped, paste0(name, "$", colnames(rd)[is_matrix]))
            rowData(se) <- rd[, !is_matrix, drop = FALSE]
        }
        if (!identical(rownames(se), key)) {
            rowData(se)$mulink_feature_id <- rownames(se)
            rownames(se) <- key
        }
        se
    }, experiments(object), names(object), keys, SIMPLIFY = FALSE)

    if (length(dropped))
        warning("Dropped matrix-valued rowData column(s) '",
                paste(dropped, collapse = "', '"), "', which MuData cannot ",
                "encode. The links they define are written to .varp.")

    ## Renaming and dropping invalidate the AssayLinks, whose hits still
    ## reference the old row names. writeH5MU() does not read them and the graph
    ## has already been extracted, so the slot is assigned directly to skip the
    ## validity check that would reject this write-only object.
    object@ExperimentList <- ExperimentList(prepared)

    MuData::writeH5MU(object, path)
    write_var_index(path, index$key)

    ## Read the index back rather than assuming it, so that .varp is aligned to
    ## what the file actually contains.
    h5 <- H5Fopen(path, flags = "H5F_ACC_RDONLY", native = FALSE)
    var_names <- .var_names(h5)
    H5Fclose(h5)

    if (!setequal(var_names, index$key))
        stop("The global .var index in '", path, "' does not match the ",
             "expected feature axis.")

    write_varp(path, feature_mapping_key,
               feature_mapping[var_names, var_names, drop = FALSE])
    invisible(path)
}
