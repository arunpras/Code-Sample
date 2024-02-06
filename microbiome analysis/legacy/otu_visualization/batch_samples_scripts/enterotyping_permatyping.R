# ENTEROTYPING & PERMATYPING
# Written By Arun Manoharan
print("enterotyping and permatyping")

# This script takes in a pyloseq object which must contain an otu table and a taxonomy table and 
# performs enterotyping and permatyping (clustering) at the genus level. Outputs are enterotype and permatype assignments,
# as well as ordination coordinates for both enterotypes and permatypes.

# test if correct number of arguments, if not return an error
NUM_ARGS <- 3
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != NUM_ARGS) {
  stop(paste("enterotyping_permatyping.R requires", NUM_ARGS, "arguments:
             PROJECT_PATH: absolute path to project directory
             PERMATYPING_NUM_PERMS: number of permutations of data to calculate enterotypes for; 1000 was used in Knomics paper
             PERMATYPING_STABILITY_CUTOFF: threshold for sample stability score; 0.56 used in Knomics paper
             
             This script requires a phyloseq object has been created from the otu and taxonomy assignment data.
             That is done by otu_visualization.R script
             phyloseq_object.RDS must be located in directory PROJECT_PATH/otu_visualization_data
             
             Example Usage:
             Rscript enterotyping_permatyping.R
             /home/arun/workdir/test_datasets/qiita_single_fastq_test
             1000
             0.56", sep=" "), call. = F)
}

# requirements for enterotyping and permatyping
library(phyloseq)
library(readr)
library(reshape2)
library(cluster)
library(clusterSim)
library(ade4)

# read command line input
PROJECT_PATH <- args[1] # absolute path to project directory
PERMATYPING_NUM_PERMS <- as.numeric(args[2]) # number of permutations of data to calculate enterotypes for; 1000 was used in Knomics paper
PERMATYPING_STABILITY_CUTOFF <- as.numeric(args[3]) # threshold for sample stability score; 0.56 used in Knomics paper

# create output directory
OTU_VIS_LOC <- paste(PROJECT_PATH, "otu_visualization_data", sep="/")
OUTPUT_LOC <- paste(PROJECT_PATH, "batch_visualization_data", "permatyping", sep="/")
dir.create(paste(PROJECT_PATH, "batch_visualization_data", sep="/"))
dir.create(OUTPUT_LOC)

# create log file to document program progress
cat(paste(Sys.time(), "Reading command line input and creating output directories completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

###########################################################################################
# load phyloseq object
cat(paste(Sys.time(), "Loading phyloseq object started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
ps_object <- readRDS(file = paste(OTU_VIS_LOC, "phyloseq_object.RDS", sep="/"))
cat(paste(Sys.time(), "Loading phyloseq object completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# get taxonomic composition
cat(paste(Sys.time(), "Getting genus level abundance from phyloseq object started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
data_tax_glom_g <- tax_glom(ps_object, 'g', NArm = F)
data_tax_glom_rarefied_g <- transform_sample_counts(data_tax_glom_g, function(OTU) OTU/sum(OTU) )

# get taxonomic composition in correct format
# get rarefied otu abundances not combined with taxonomic information - found in data_tax_glom_rarefied_g (created during taxonomic distribution step above)
# genus level enterotyping is standard (see: http://enterotyping.embl.de/)
data_enterotyping <- otu_table(data_tax_glom_rarefied_g)
cat(paste(Sys.time(), "Getting genus level abundance from phyloseq object completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# define enterotyping and permatyping functions
cat(paste(Sys.time(), "Defining enterotyping and permatyping functions started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
# define distance metric functions as given in tutorial: http://enterotype.embl.de/enterotypes.html
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}

# define clustering function as given in tutorial: http://enterotype.embl.de/enterotypes.html
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

# define enterotyping function for permatyping
get_enterotypes = function(ent_data, is_get_orig_enterotype) {
  # calculate distance between samples
  ent_data.dist=dist.JSD(ent_data)
  # if getting original enterotypes, determine optimal number of clusters (choose number of clusters with highest CH index value)
  if (is_get_orig_enterotype) {
    nclusters=NULL
    for (k in 1:min(20, length(colnames(ent_data)) - 1)) {
      if (k==1) {
        nclusters[k]=NA
      } else {
        ent_data.cluster_temp=pam.clustering(ent_data.dist, k)
        nclusters[k]=index.G1(t(ent_data),ent_data.cluster_temp,  d = ent_data.dist, centrotypes = "medoids")
      }
    }
    num_clusters = which.max(nclusters)
    # otherwise, use pre-determined optimal number of clusters
  } else {
    num_clusters = orig_num_clusters
  }
  # print(paste("number of clusters: ", num_clusters, sep=""))
  # cluster samples
  ent_data.cluster=pam.clustering(ent_data.dist, k=num_clusters)
  nclusters = index.G1(t(ent_data), ent_data.cluster, d = ent_data.dist, centrotypes = "medoids")
  return(list("clusters" = ent_data.cluster, "n" = num_clusters, "dists" = ent_data.dist))
}
cat(paste(Sys.time(), "Defining enterotyping and permatyping functions completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# get original enterotypes
cat(paste(Sys.time(), "Getting original enterotypes started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
results_list = get_enterotypes(data_enterotyping, is_get_orig_enterotype=T)
orig_ent_data.cluster = results_list$clusters
orig_num_clusters = results_list$n
orig_ent_data.dist = results_list$dists

orig_ET = as.data.frame(cbind(as.character(colnames(data_enterotyping)), orig_ent_data.cluster))
ET_totals = as.data.frame(table(orig_ent_data.cluster))
cat(paste(Sys.time(), "Getting original enterotypes completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# write out original enterotypes
cat(paste(Sys.time(), "Writing out original enterotype data to file started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
# rename columns from "V1" and "orig_ent_data.cluster" to more informative names
colnames(orig_ET) <- c("sample_name", "cluster")
write.table(orig_ET, file=paste(OUTPUT_LOC, 'original_enterotypes.txt', sep="/"), sep="\t", row.names = F, quote = F)

# do ordination with original enterotypes
ent.pcoa=dudi.pco(orig_ent_data.dist, scannf=F, nf=2)
# write out original enterotype PCOA coordinates
# rename axes from A1 and A2 to more descriptive names
ent_coordinates <- ent.pcoa$li
colnames(ent_coordinates) <- c("PCoA Axis 1", "PCoA Axis 2")
write.table(ent_coordinates, file=paste(OUTPUT_LOC, 'orig_ent_coordinates.txt', sep="/"), sep="\t", row.names = F, quote = F)
cat(paste(Sys.time(), "Writing out original enterotype data to file completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# do permatyping
cat(paste(Sys.time(), "Permutations of enterotyping started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
permatype_clusters = matrix(, nrow = dim(data_enterotyping)[[2]], ncol = PERMATYPING_NUM_PERMS)
all_samples = as.data.frame(colnames(data_enterotyping))

for (i in seq(1, PERMATYPING_NUM_PERMS)) {
  # subsample 50% of samples from each original enterotype
  bootstrap_samples = vector()
  # for every enterotype cluster, sample 50% of samples without replacement
  for (j in seq(1, orig_num_clusters)) {
    # get indices of samples in enterotype j (character type because integer gets misinterpreted by sample function)
    enterotype_j_indices = c(rownames(orig_ET)[which(orig_ET$cluster == j)])
    # sample 50% of the samples in enterotype j, record those indices
    sampled_indices = sample(enterotype_j_indices, size = as.integer(ET_totals[j, 'Freq']/2) + 1)
    bootstrap_samples = c(bootstrap_samples, as.integer(sampled_indices))
  }
  ent_data = data_enterotyping[, bootstrap_samples]
  # re-rarify subsample table
  ent_data = sweep(ent_data, 2, colSums(ent_data), `/`)
  
  # repeat enterotyping with these samples
  results_list = get_enterotypes(ent_data, is_get_orig_enterotype=F)
  new_clusters = results_list$clusters
  orig_num_clusters = results_list$n
  ent_data.dist = results_list$dists
  
  new_ET = as.data.frame(cbind(as.character(colnames(ent_data)), new_clusters))
  all_samples_ET = merge(x = all_samples, y = new_ET, by.x = 'colnames(data_enterotyping)', by.y = 'V1', all.x = T)
  
  permatype_clusters[,i] = all_samples_ET$new_clusters
}

perm_clusters = as.data.frame(permatype_clusters)
rownames(perm_clusters) = colnames(data_enterotyping)
cat(paste(Sys.time(), "Permutations of enterotyping completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# build permatyping distance matrix
cat(paste(Sys.time(), "Building sample by sample distance matrix based on permutations of enterotyping started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
ent_dists = matrix(, nrow = dim(data_enterotyping)[[2]], ncol = dim(data_enterotyping)[[2]])
for (i in seq(1, dim(perm_clusters)[[1]])) {
  for (j in seq(1, dim(perm_clusters)[[1]])) {
    comparison_table = table(perm_clusters[i, ] == perm_clusters[j, ])
    Dij = comparison_table['TRUE']/ sum(comparison_table)
    ent_dists[i, j] = Dij
  }
}

ent_dists = as.data.frame(ent_dists)
rownames(ent_dists) = colnames(data_enterotyping)
colnames(ent_dists) = colnames(data_enterotyping)
cat(paste(Sys.time(), "Building sample by sample distance matrix based on permutations of enterotyping completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

# determine which samples are stable according to permatyping
cat(paste(Sys.time(), "Determining stable samples based on sample distance to enterotype medoids started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
# find original cluster mediods given dissimilarity matrix (orig_ent_data.dist)
medoids = as.vector(pam(orig_ent_data.dist, orig_num_clusters, diss=TRUE))$medoids
print(paste("medoids: ", medoids, sep=""))

# keep only samples with stability scores > PERMATYPING_STABILITY_CUTOFF
# stability score defined in correspondance with author as max(c(dist_ent1, ..., dist_entN)) / sum(c(dist_ent1, ..., dist_entN))
colMax <- function(data) sapply(data, max, na.rm = T)
stability_scores = colMax(ent_dists[medoids, ])/colSums(ent_dists[medoids, ], na.rm = T)
stable_samples = names(stability_scores[which(stability_scores > PERMATYPING_STABILITY_CUTOFF)])
cat(paste(Sys.time(), "Determining stable samples based on sample distance to enterotype medoids completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

print(paste("# stable samples: ", length(stable_samples), sep=""))
print(paste("# enterotypes: ", orig_num_clusters, sep=""))

# generate enterotypes using only stable samples
cat(paste(Sys.time(), "Generating enterotypes using only stable samples (now called permatypes) started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
print("getting permatypes with stable samples only")
ent_data = data_enterotyping[, which(colnames(data_enterotyping) %in% stable_samples)]
cat(paste(Sys.time(), "Generating enterotypes using only stable samples (now called permatypes) completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")

cat(paste(Sys.time(), "Checking if number of stable samples > number of clusters, writing out permatyping results started.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")
if (length(stable_samples) <= orig_num_clusters) {
  print("number of stable samples <= number of enterotypes, permatyping fails")
} else {
  results_list = get_enterotypes(ent_data, is_get_orig_enterotype=F)
  final_ent_data.cluster = results_list$clusters
  num_clusters = results_list$n
  final_ent_data.dist = results_list$dists
  
  final_ET = as.data.frame(cbind(as.character(colnames(ent_data)), final_ent_data.cluster))
  final_ET_totals = as.data.frame(table(final_ent_data.cluster))
  
  # do ordination with final permatypes
  final_ent.pcoa=dudi.pco(final_ent_data.dist, scannf=F, nf=2)
  
  # write out final permatypes
  # rename columns from "V1" and "final_ent_data.cluster" to more informative names
  colnames(final_ET) <- c("sample_name", "cluster")
  write.table(final_ET, file=paste(OUTPUT_LOC, 'final_permatypes.txt', sep="/"), sep="\t", row.names = F, quote = F)
  # write out final permatype PCOA coordinates
  # rename axes from A1 and A2 to more descriptive names
  final_ent_coordinates <- final_ent.pcoa$li
  colnames(final_ent_coordinates) <- c("PCoA Axis 1", "PCoA Axis 2")
  write.table(final_ent_coordinates, file=paste(OUTPUT_LOC, 'final_perm_coordinates.txt', sep="/"), sep="\t", row.names = F, quote = F)
}
cat(paste(Sys.time(), "Checking if number of stable samples > number of clusters, writing out permatyping results completed.\n", sep=" "), file=paste(OUTPUT_LOC, "enterotyping_permatyping_log.txt", sep="/"), append=T, sep="/n")