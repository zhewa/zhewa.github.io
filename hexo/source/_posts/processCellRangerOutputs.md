---
title: An R script to process Cell Ranger output files
date: 2019-05-01 10:29:03
tags:
- Single Cell
- Script
- Cell Ranger
- R
---
First install the required packages `data.table`, `Matrix`, and `SingleCellExperiment`.
{% codeblock lang:objc %}
if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
    install.packages("Matrix")
}

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("SingleCellExperiment")
}
{% endcodeblock %}<!-- more -->The definition of function `processCellRangerOutput` is shown below. This function reads the count files (in sparse matrix format), the cell barcodes, and the feature annotations for each sample, writes each of the count matrices as a tab separated text file, then combines the count matrices as one for all samples, and finally writes it as a tab separated text file as well as a `SingleCellExperiment` object to specified path.

Note that this process requires a lot of memory depending on the size of the data. It might fail if the free memory on the computer is not enough for allocation or the data is too large to fit in the RAM.
{% codeblock processCellRangerOutput.R lang:r https://github.com/zhewa/usefulscripts/blob/master/processCellRangerOutput.R %}
#' @title Process Cell Ranger output
#' @description Read the filtered matrices for all samples from Cell Ranger
#'  output. Combine them into one big expression matrix and write it to
#'  \code{outDir}. A \link[SingleCellExperiment]{SingleCellExperiment} object is
#'  also saved as a `.RDS` file.
#' @param cellRangerDir The root path where Cell Ranger was run. This folder
#'  should contain sample specific folders.
#' @param samples A vector of sample names. Must be the same as the folder names
#'  of the samples.
#' @param project The project name. Can be any string. It is
#'  added to the filenames of the combined count matrix and the
#'  \code{SingleCellExperiment} object.
#' @param writeIndivisualSampleCounts Boolean. Whether to write sample specific
#'  count matrices to output files or not. Default \code{TRUE}.
#' @param writeAllSampleCounts Boolean. Whether to write combined count matrix
#'  to output file or not. Default \code{TRUE}.
#' @param writeSCE Boolean. Whether to save the \code{SingleCellExperiment}
#'  object containing the combined count matrix to output file or not. Default
#'  \code{TRUE}.
#' @param sep The field separator string. Values within each row are separated
#'  by this string when writing output files.
#' @param outDir The path to output location. Default to current working
#'  directory.
#' @param outPrefix Prefix to output files. Default to current date.
#' @param overwrite Boolean. Whether to overwrite output files or not. Default
#'  \code{FALSE}.
#' @param cellRangerOuts The intermidiate path to filtered feature count files
#'  saved in sparce matrix format. Reference genome names might need to be
#'  appended if reads were mapped to multiple genomes when running Cell Ranger
#'  pipeline. Default \code{"outs/filtered_feature_bc_matrix/"}.
#' @param gzipped Boolean. Whether the Cell Ranger output files were gzip
#'  compressed or not. This is true subsequent to Cell Ranger 3.0. Default
#'  \code{TRUE}.
#' @return A \code{SingleCellExperiment} object containing the combined count
#'  matrix, the feature annotations, and the cell annotation.
processCellRangerOutput <- function(
    cellRangerDir,
    samples,
    project = strtrim(paste0(samples, collapse = ""), 5),
    writeIndivisualSampleCounts = TRUE,
    writeAllSampleCounts = TRUE,
    writeSCE = TRUE,
    sep = "\t",
    outDir = ".",
    outPrefix = format(Sys.time(), "%Y%m%d"),
    overwrite = FALSE,
    cellRangerOuts = "outs/filtered_feature_bc_matrix/",
    gzipped = TRUE) {

    fileNames <- c(matrix = "matrix.mtx",
        features = "features.tsv",
        barcodes = "barcodes.tsv")
    if (isTRUE(gzipped)) {
        fileNames <- vapply(fileNames,
            FUN = function(i) {
                paste0(i, ".gz")
            },
            FUN.VALUE = character(1))
    }

    res <- vector("list", length = length(samples))
    coldata <- vector("list", length = length(samples))

    for (i in seq_along(samples)) {
        mm <- Matrix::readMM(gzfile(file.path(cellRangerDir,
            samples[i], cellRangerOuts, fileNames[["matrix"]])))

        fe <- read.table(gzfile(file.path(cellRangerDir,
            samples[i], cellRangerOuts, fileNames[["features"]])), sep = "\t")

        bc <- read.table(gzfile(file.path(cellRangerDir,
            samples[i], cellRangerOuts, fileNames[["barcodes"]])), sep = "\t")
        coldata[[i]] <- rep(samples[i], nrow(bc))

        ma <- as.matrix(mm)
        colnames(ma) <- paste(bc[[1]], samples[i], sep = "_")
        rownames(ma) <- fe[[1]]
        res[[i]] <- ma

        if (isTRUE(writeIndivisualSampleCounts)) {
            fn <- file.path(outDir,
                paste0(outPrefix, "_", samples[i], ".txt"))
            if (isFALSE(overwrite) & file.exists(fn)) {
                warning("File already exists! Skip writting file ", fn)
            } else {
                data.table::fwrite(ma,
                    file = fn,
                    sep = sep,
                    row.names = TRUE)
            }
        }
    }

    expr <- do.call(cbind, res)
    exprdt <- data.table::as.data.table(expr, keep.rownames = TRUE)
    colnames(exprdt)[1] <- "geneID"

    if (isTRUE(writeAllSampleCounts)) {
        fn <- file.path(outDir,
            paste0(outPrefix, "_", project, "_allcounts.txt"))
        if (isFALSE(overwrite) & file.exists(fn)) {
            warning("File already exists! Skip writting file ", fn)
        } else {
            data.table::fwrite(exprdt,
                file = fn,
                sep = sep,
                row.names = TRUE)
        }
    }

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = expr))
    SummarizedExperiment::rowData(sce) <- fe
    colnames(SummarizedExperiment::rowData(sce)) <- c("geneID",
        "genename", "genebiotype")
    coldata <- unlist(coldata)
    SummarizedExperiment::colData(sce)["sample"] <- data.frame(sample = coldata)

    if (isTRUE(writeSCE)) {
        fn <- file.path(outDir,
            paste0(outPrefix, "_", project, "_allcounts_sce.RDS"))
        if (isFALSE(overwrite) & file.exists(fn)) {
            warning("File already exists! Skip writting file ", fn)
        } else {
            saveRDS(sce, file = fn)
        }
    }

    return(sce)
}

{% endcodeblock %}