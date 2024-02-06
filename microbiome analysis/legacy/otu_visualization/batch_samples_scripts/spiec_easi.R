# SPIEC-EASI: TAXA CO-OCCURANCE
# Written By Arun Manoharan
print("analysing taxa co-occurance")

# This script takes in a phyloseq object which must contain an otu table and returns a cooccurrance table
# with whether or not 2 taxa co-occur across samples, as well as an abundace table to scale nodes of the 
# cooccurance graph by if desired. 

# test if correct number of arguments, if not return an error
NUM_ARGS <- 1
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != NUM_ARGS) {
  stop(paste("spiec_easi.R requires", NUM_ARGS, "arguments:
             PROJECT_PATH: absolute path to project directory

             This script requires a phyloseq object has been created from the otu and taxonomy assignment data.
             That is done by otu_visualization.R script
             phyloseq_object.RDS must be located in directory PROJECT_PATH/otu_visualization_data
             
             Example Usage:
             Rscript spiec_easi.R
             /home/arun/workdir/test_datasets/qiita_single_fastq_test", sep=" "), call. = F)
}

# requirements for cooccurrance networks - spiec easi
library(phyloseq)
library(readr)
library(reshape2)
library(SpiecEasi)

# spieceasi:
#   manual: https://stamps.mbl.edu/images/e/e9/SpiecEasi_tutorial_08102016.pdf
#   installation: 
#     install.packages("devtools")
#     devtools::install_github("zdk123/SpiecEasi", force=TRUE)

# read command line input
PROJECT_PATH <- args[1] # absolute path to project directory

# create output directory
OTU_VIS_LOC <- paste(PROJECT_PATH, "otu_visualization_data", sep="/")
OUTPUT_LOC <- paste(PROJECT_PATH, "batch_visualization_data", "cooccurance_spiec_easi", sep="/")
dir.create(paste(PROJECT_PATH, "batch_visualization_data", sep="/"))
dir.create(OUTPUT_LOC)

# create log file to document program progress
cat(paste(Sys.time(), "Reading command line input and creating output directories completed.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")

###########################################################################################
# load phyloseq object
cat(paste(Sys.time(), "Loading phyloseq object started.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")
ps_object <- readRDS(file = paste(OTU_VIS_LOC, "phyloseq_object.RDS", sep="/"))
cat(paste(Sys.time(), "Loading phyloseq object comleted.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")

# display original otu table dimensions
cat(paste(Sys.time(), "Filtering otu table started.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")
print("original otu table dimensions (taxa by samples):")
print(dim(otu_table(ps_object)))

# determine first quartile read depth across samples
first_quartile_read_depth = quantile(sample_sums(ps_object))[2]

# filter phyloseq object to remove OTUS present in fewer than 37% of samples (threshold used in spiec-easi paper)
data_se_filtered = filter_taxa(ps_object, function(x) sum(x > 0) > (0.37*length(x)), TRUE)

# filter phyloseq object to remove samples if sequencing depth fell below 1st quartile (threshold used in spiec-easi paper)
data_se_filtered = prune_samples(sample_sums(data_se_filtered) > first_quartile_read_depth, data_se_filtered)

# display final otu table dimensions
print("otu table dimensions (taxa by samples) after filtering samples with read depth < first quartile depth and taxa present in fewer than 37% of samples:")
print(dim(otu_table(data_se_filtered)))
cat(paste(Sys.time(), "Filtering otu table completed.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")

# run speic-easi to generate adjacency matrix
cat(paste(Sys.time(), "Calling spiec.easi to generate adjacency matrix started.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")
se.mb <- spiec.easi(data_sse_filtered, method='mb', lambda.min.ratio=1e-2, nlambda=20, icov.select.params=list(rep.num=50))
cat(paste(Sys.time(), "Calling spiec.easi to generate adjacency matrix completed.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")

# output co-occurance (adjacent nodes) table
cat(paste(Sys.time(), "Writing out co-occurance table started.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")
adjacency_table <- as.data.frame(as.matrix(se.mb$refit))
colnames(adjacency_table) <- taxa_names(data_se_filtered)
rownames(adjacency_table) <- taxa_names(data_se_filtered)
write.table(adjacency_table, file = paste(OUTPUT_LOC, "coocurrance_table.txt", sep="/"), sep="\t", quote = F, row.names = T, col.names = NA)
cat(paste(Sys.time(), "Writing out co-occurance table completed.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")

# output abundance table
cat(paste(Sys.time(), "Writing out abundance table started.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")
write.table(otu_table(data_se_filtered), file=paste(OUTPUT_LOC, "abundance_matrix.txt", sep="/"), sep="\t", quote = F, row.names = T, col.names = NA)
cat(paste(Sys.time(), "Writing out abundance table completed.\n", sep=" "), file=paste(OUTPUT_LOC, "speic_easi_log.txt", sep="/"), append=T, sep="/n")
