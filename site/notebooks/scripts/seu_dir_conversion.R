# Function to export Seurat object to directory for AnnData loading.
# Excludes graphs (neighbor graphs).
# Compatible with Seurat v4 and v5 (checks assay class/slots).
# Assumes 'RNA' assay; adjust if needed.

# R >= 4.1; Seurat 3.x/4.x/5.x compatible
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

.as_dgC <- function(m) {
    if (inherits(m, "dgCMatrix")) return(m)
    if (is(m, "sparseMatrix")) return(as(m, "dgCMatrix"))
    Matrix(m, sparse = TRUE)
}

.collect_assay_layers <- function(a) {
    # Return a named list of matrices for export; supports Seurat 3/4 and 5+
    out <- list()
    # Seurat 5+: layers slot
    if ("layers" %in% slotNames(a)) {
        lyr <- a@layers
        if (!is.null(lyr) && length(lyr) > 0) {
            for (nm in names(lyr)) out[[nm]] <- lyr[[nm]]
        }
    }
    # Back-compat canonical slots (present in 3/4 and still valid in 5)
    if ("counts" %in% slotNames(a) && !is.null(a@counts))       out[["counts"]] <- a@counts
    if ("data" %in% slotNames(a)   && !is.null(a@data))         out[["data"]]   <- a@data
    out
}

export_seurat_to_anndata_dir <- function(
        seu_obj,
        output_parent_dir,
        sample,
        assay = DefaultAssay(seu_obj)
) {
    # --- 0. Paths ---
    sample_dir <- file.path(output_parent_dir, sample)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    
    # --- 1. Basic metadata ---
    stopifnot(assay %in% names(seu_obj@assays))
    a <- seu_obj@assays[[assay]]
    
    # Cells (obs)
    pmeta <- seu_obj@meta.data
    if (is.null(rownames(pmeta))) stop("meta.data must have rownames (cell barcodes).")
    
    # Features (var) from the chosen assay
    feat_names <- rownames(a)
    fmeta <- data.frame(gene_names = feat_names, row.names = feat_names, check.names = FALSE)
    
    # --- 2. Layers: version-agnostic collection ---
    layers_list <- .collect_assay_layers(a)
    if (length(layers_list) == 0) warning("No layers found in assay '", assay, "'. Nothing to write for matrices.")
    
    for (nm in names(layers_list)) {
        mat <- layers_list[[nm]]
        # Ensure dimnames exist (genes x cells)
        if (is.null(rownames(mat))) rownames(mat) <- feat_names
        if (is.null(colnames(mat))) colnames(mat) <- colnames(seu_obj)
        Matrix::writeMM(.as_dgC(mat), file = file.path(sample_dir, paste0(nm, ".mtx")))
    }
    
    # --- 3. Write obs / var ---
    write.table(
        pmeta, file = file.path(sample_dir, "obs.tsv"),
        sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
    )
    write.table(
        fmeta, file = file.path(sample_dir, "var.tsv"),
        sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
    )
    
    # --- 4. Reductions: export all available embeddings ---
    reds <- names(seu_obj@reductions)
    for (r in reds) {
        emb <- tryCatch(Embeddings(seu_obj, reduction = r), error = function(e) NULL)
        if (!is.null(emb) && nrow(emb) > 0) {
            emb_df <- as.data.frame(emb, check.names = FALSE)
            # Save as cells x dims; rownames = cell IDs
            write.table(
                emb_df,
                file = file.path(sample_dir, paste0(r, "_embeddings.tsv")),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
            )
        }
    }
    
    invisible(sample_dir)
}


# Example usage:
# library(Seurat)
# source("seu_dir_conversion.R")
# data(pbmc_small)  # Loads the Seurat object (v4 compatible)
# export_seurat_to_anndata_dir(pbmc_small, "./data", "pbmc_small_test")





# Read a TSV written by pandas (index in the first column)
.read_tsv <- function(path) {
    df <- tryCatch(
        as.data.frame(fread(path, sep = "\t", header = TRUE), stringsAsFactors = FALSE),
        error = function(e) read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    )
    rownames(df) <- df[[1]]
    df[[1]] <- NULL
    df
}

# Read MatrixMarket (genes x cells) and return dgCMatrix
.read_mtx_gxc <- function(path) {
    m <- Matrix::readMM(path)                # returns dgTMatrix
    m <- as(m, "CsparseMatrix")              # preferred, avoids deprecation warning
    storage.mode(m@x) <- "double"            # ensure numeric
    m
}

.set_data_layer <- function(seu, assay, mat) {
    so_ver <- tryCatch(utils::packageVersion("SeuratObject"), error = function(e) "0.0.0")
    if (!inherits(so_ver, "package_version")) so_ver <- "0.0.0"
    if (so_ver >= "5.0.0") {
        # New API: use layers
        seu <- SetAssayData(seu, assay = assay, layer = "data", new.data = mat)
    } else {
        # Old API: data slot
        seu <- SetAssayData(seu, assay = assay, slot  = "data", new.data = mat)
    }
    seu
}

.reduction_key <- function(nm_clean) {
    base <- gsub("[^A-Za-z0-9]", "", toupper(nm_clean))  # keep only A–Z,0–9
    if (nchar(base) == 0) base <- "DR"
    paste0(base, "_")
}


# Main loader
import_seurat_from_dir <- function(dir, assay = "RNA") {
    dir <- normalizePath(dir, mustWork = TRUE)
    
    # --- obs / var ---
    obs <- .read_tsv(file.path(dir, "obs.tsv"))   # cells as rows
    var <- .read_tsv(file.path(dir, "var.tsv"))   # genes as rows
    
    # --- matrices (genes x cells on disk) ---
    counts_fp <- file.path(dir, "counts.mtx")
    data_fp   <- file.path(dir, "data.mtx")
    has_counts <- file.exists(counts_fp)
    has_data   <- file.exists(data_fp)
    if (!has_counts && !has_data) {
        stop("No counts.mtx or data.mtx found in: ", dir)
    }
    
    # counts: prefer counts.mtx, else fallback to data.mtx
    counts_gxc <- .read_mtx_gxc(if (has_counts) counts_fp else data_fp)
    rownames(counts_gxc) <- rownames(var)
    colnames(counts_gxc) <- rownames(obs)
    
    seu <- CreateSeuratObject(counts = counts_gxc, assay = assay, meta.data = obs)
    
    # normalized data slot, if present
    if (has_data) {
        data_gxc <- .read_mtx_gxc(data_fp)
        rownames(data_gxc) <- rownames(var)
        colnames(data_gxc) <- rownames(obs)
        seu <- .set_data_layer(seu, assay = assay, mat = data_gxc)
    }
    
    # --- ALL embeddings: files named like 'X_<name>_embeddings.tsv' (or '<name>_embeddings.tsv') ---
    emb_files <- list.files(dir, pattern = "_embeddings\\.tsv$", full.names = TRUE)
    for (ef in emb_files) {
        base <- sub("_embeddings\\.tsv$", "", basename(ef))  # e.g., 'X_pca', 'X_umap', 'lsi'
        nm_clean <- sub("^X_", "", base)                      # strip leading 'X_' if present
        emb <- .read_tsv(ef)                                  # cells x dims
        # align to Seurat cell order
        emb <- emb[colnames(seu), , drop = FALSE]
        key <- .reduction_key(nm_clean)
        dr  <- CreateDimReducObject(embeddings = as.matrix(emb), key = key, assay = assay)
        seu[[nm_clean]] <- dr
    }
    
    seu
}

# Example usage:
# library(Seurat)
# source("seu_dir_conversion.R")
# seu_res = import_seurat_from_dir("./data/pbmc_small_test_wt_concord/")



