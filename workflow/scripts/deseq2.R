# activate the conda env 'deseq2' before running this script and open R terminal

# DESeq2 script to create and MDS and heatmap.



# https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/DESeq2.html

library(readr)
library(DESeq2)
library(tximport)
library(apeglm)

### set path to directory containing the salmon output per sample

dir <- 'rnaseq_AS_Esup/star_salmon'

### set path to output directory and create if it doesn't exist

output_dir <- 'deseq2_results_AS_Esup'
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

### get sample names from directory

samples <- gsub("/quant.sf", "", list.files(dir, recursive=TRUE, pattern="quant.sf"))

### load metadata

# read all metadata
md <- read.table("additional_data/sample_list_metadata.txt", sep="\t", header=TRUE)

# create a column for sample name i.e. group_replicate
md$samples <- paste(md$group, md$replicate, sep="_")

# check metadata
head(md)

### create tx2gene mapping

tx2gene <- read.table('/home/owhite/my_projects/data_krill_chromosome_annotations/lx.sliding_filtered.bed',
                   sep="\t", header=FALSE, col.names=c("chr","start","end","gene","score","strand"))
tx2gene$GENEID <- tx2gene$gene
tx2gene$TXNAME <- sub("_exon\\d+$", "", tx2gene$gene)
tx2gene <- tx2gene[,c("TXNAME","GENEID")]
head(tx2gene)

### define groups to compare

control <- "AS_Esup_0C_0hr"
treatments <- c("AS_Esup_3C_3hr", "AS_Esup_3C_6hr", "AS_Esup_6C_3hr","AS_Esup_6C_6hr")

### loop control treatment combinations
for (treatment in treatments) {

    # manually set control and treatment for testing
    # treatment = "AS_Esup_3C_3hr"

    cat(paste("Comparing", control, "vs", treatment, "\n"))

    # set the names of the groups to be compared
    group_selected <- c(control, treatment)

    # filter the metadata to only include the groups of interest
    cat("Filtering metadata for selected groups\n")
    md_filtered <- md[ md$group %in% group_selected, ]

    # subset md to only include samples that exist in the quant.sf files
    md_filtered <- md_filtered[ md_filtered$samples %in% samples, ]

    # set group as factor with levels in desired order
    md_filtered$group <- factor(md_filtered$group, levels=group_selected)

    # show some summary info
    cat(paste("   Total samples", length(md_filtered$samples), "\n"))
    cat(paste("   Total groups", length(unique(md_filtered$group)), "\n"))
    cat(paste("   Control samples", sum(md_filtered$group == control), "\n"))
    cat(paste("   Treatment samples", sum(md_filtered$group == treatment), "\n"))

    # get file paths to quant.sf files
    files <- file.path(dir, md_filtered$samples, "quant.sf")

    # import transcript-level quantifications and summarize to gene-level
    cat("Importing transcript-level quantifications and summarizing to gene-level\n")
    txi <- tximport(files, type="salmon", tx2gene=tx2gene)

    # create DESeq2 dataset
    cat("Creating DESeq2 dataset\n")
    dds <- DESeqDataSetFromTximport(txi, colData=md_filtered, design= ~ group)

    nrow(dds)

    # check for genes with no counts across all samples
    cat(paste(sum(rowSums(counts(dds)) == 0), "of", nrow(dds), "genes have no coverage:", "\n"))
    # get the smallest group size
    smallestGroupSize <- min(table(md_filtered$group))
    # keep genes with at least 10 counts in at least smallestGroupSize samples
    keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    cat("Prefiltering genes, keeping genes with at least 10 counts in at least", smallestGroupSize, "samples\n")
    cat(paste("Keeping ", sum(keep), "out of", nrow(dds), "genes"), "\n")
    # filter dds
    dds <- dds[keep,]

    # run DESeq2
    cat("Running DESeq2\n")
    dds <- DESeq(dds)
    # get results with default alpha = 0.05
    res <- results(dds, alpha=0.05)

    # shrink log2 fold changes for more accurate estimates
    cat("Applying LFC shrinkage\n")
    resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

    # show summary of results
    cat("DESeq2 results summary:\n")
    summary(res)

    # how many genes are differentially expressed at FDR < 0.05?
    cat("   Number of DE genes at FDR < 0.05:", sum(res$padj < 0.05, na.rm=TRUE), "\n")

    # sort results by adjusted p-value
    res <- res[order(res$padj), ]

    # write results to file
    cat("Writing results to file\n")
    output_file <- file.path(output_dir, paste0("DESeq2_results_", control, "_vs_", treatment, ".tsv"))
    write.table(as.data.frame(res), file=output_file, sep="\t", quote=FALSE, row.names=TRUE)

    # print MA plot
    output_file_ma_plot <- file.path(output_dir,  paste0("DESeq2_results_MA_plot_", control, "_vs_", treatment, ".png"))
    png(output_file_ma_plot, res=100, height=6, width=8, units="in")
    plotMA(res, ylim=c(-2,2))
    dev.off()

    # print MA plot with LFC shrinkage
    output_file_ma_plot_LFC <- file.path(output_dir,  paste0("DESeq2_results_MA_plot_LFCshrunken_", control, "_vs_", treatment, ".png"))
    png(output_file_ma_plot_LFC, res=100, height=6, width=8, units="in")
    plotMA(resLFC, ylim=c(-2,2))
    dev.off()

}

    # # get normalized counts for DE genes only
    # degs <- rownames(res)[ which(res$padj < 0.05) ]
    #
    # # Only proceed if we have DE genes
    # if(length(degs) > 0) {
    #     cat(paste("Extracting normalized counts for", length(degs), "DE genes\n"))
    #     counts_deg <- as.data.frame(counts(dds[degs, ], normalized=TRUE))
    #
    #     # set column names to sample names
    #     colnames(counts_deg) <- md_filtered$samples
    #     output_file_degs_counts <- file.path(output_dir, paste0("DESeq2_DEG_counts_", control, "_vs_", treatment, ".tsv"))
    #     write.table(counts_deg, file=output_file_degs_counts, sep="\t", quote=FALSE, row.names=TRUE)
    # } else {
    #     cat("No significantly differentially expressed genes found (padj < 0.05)\n")
    # }
