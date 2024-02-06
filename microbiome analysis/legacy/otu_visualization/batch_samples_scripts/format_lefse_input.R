###########################################################################################
# Written by Arun Manoharan (arun@primediscoveries.com)
# This script formats metadata for lefse biomarker analysis

# https://twbattaglia.gitbooks.io/introduction-to-qiime/content/lefse.html
# https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial
###########################################################################################

# take in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# test if correct number of arguments, if not return an error
NUM_ARGS <- 5
if (length(args) != NUM_ARGS) {
  stop(paste("format_lefse_input.R requires", NUM_ARGS, "arguments: 
       PROJECT_PATH: absolute path to project directory
       METADATA_LOC: relative path from project directory to metadata file, including file name and extension
       LEFSE_SUBJECT: metadata column name representing sample or patient ID ### right now must correpond to sample names in OTU table
       LEFSE_CLASS: metadata column name representing factor to determine differential abundance based on
       LEFSE_SUBCLASS: [Optional] stratefying factor under LEFSE_CLASS; Options are [<factor column name>, NA] 
       
       This script requires the phyloseq object produced by otu_visualization.R located at PROJECT_PATH/otu_visualization_data/phyloseq_object.RDS
       and a metadata file with the columns indicated in the argument descriptions above.

       Output is lefse_feature_table.txt in PROJECT_PATH/batch_visualization_data/
       
       Example Usage:
          Rscript format_lefse_input.R 
          /home/ubuntu/microbiome_downloads/soil/805
          /metadata/805_20180418-110618.txt
          sample_name
          year
          pH", sep=" "), call. = F)
}

# assign variables based on command line arguments
PROJECT_PATH <- args[1]
METADATA_LOC <- args[2]
LEFSE_SUBJECT <- args[3]
LEFSE_CLASS <- args[4]
LEFSE_SUBCLASS <- args[5]

# set filepath variables
PHYLOSEQ_PATH <- paste(PROJECT_PATH, "otu_visualization_data", "phyloseq_object.RDS", sep="/")
METATDATA_PATH <- paste(PROJECT_PATH, METADATA_LOC, sep="/")
# these directories are made by biomarker_lefse_script.sh
BATCH_DATA_DIR <- paste(PROJECT_PATH, "batch_visualization_data", sep="/")
LEFSE_OUTPUT_PATH <- paste(BATCH_DATA_DIR, "biomarker_lefse", sep="/")

# continue to write out log file
cat(paste(Sys.time(), "In format_lefse_input.R: Setting variable values completed.\n", sep=" "), file=paste(BATCH_DATA_DIR, "biomarker_lefse_script_log.txt", sep="/"), append=T, sep="/n")

# requirements
library(phyloseq)
library(readr)

# arrange otu table with taxonomic assignment as 1st column
cat(paste(Sys.time(), "In format_lefse_input.R: Adding metadata to ASV table started.\n", sep=" "), file=paste(BATCH_DATA_DIR, "biomarker_lefse_script_log.txt", sep="/"), append=T, sep="/n")

# load metadata file
metadata <- read.csv(METATDATA_PATH, sep='\t', header = T)

# load phyloseq object with feature table and taxonomic assignments
ps_object <- readRDS(PHYLOSEQ_PATH)

# get relative abundance otu table
ps_object_rel_abundance <- transform_sample_counts(ps_object, function(x) x / sum(x))
otu_tab <- otu_table(ps_object_rel_abundance)

# add user-specified metadata columns for LEFSE_CLASS and LEFSE_SUBCLASS to metadata table
# make LEFSE_CLASS column of metadata character entries
metadata[,LEFSE_CLASS] <- as.character(metadata[,LEFSE_CLASS])
# select LEFSE_CLASS entries in the order that the OTU table samples are in
lefse_class_vec <- metadata[match(colnames(otu_tab), metadata[,LEFSE_SUBJECT]), colnames(metadata)[which(colnames(metadata) == LEFSE_CLASS)]]
# add the LEFSE_CLASS vector to the top of the otu table
lefse_df <- rbind.data.frame(lefse_class_vec, otu_tab)
rownames(lefse_df)[1] <- LEFSE_CLASS

# if stratefying LEFSE_SUBCLASS is specified, do the same for subclass factor
if (LEFSE_SUBCLASS != "NA") {
  metadata[,LEFSE_SUBCLASS] <- as.character(metadata[,LEFSE_SUBCLASS])
  lefse_subclass_vec <- metadata[match(colnames(otu_tab), metadata[,LEFSE_SUBJECT]), colnames(metadata)[which(colnames(metadata) == LEFSE_SUBCLASS)]]
  lefse_df <- rbind.data.frame(lefse_subclass_vec, lefse_df)
  rownames(lefse_df)[1] <- LEFSE_SUBCLASS
}

# make column name for Feature IDs "Feature_ID" rather than empty
lefse_df <- data.frame("Feature_ID"=rownames(lefse_df), lefse_df)

# write metadata-appended otu table to file for lefse
write.table(lefse_df, file=paste(LEFSE_OUTPUT_PATH, "lefse_input_table.txt", sep="/"), sep="\t", quote = F, row.names = F)
cat(paste(Sys.time(), "In format_lefse_input.R: Adding metadata to ASV table completed.\n", sep=" "), file=paste(BATCH_DATA_DIR, "biomarker_lefse_script_log.txt", sep="/"), append=T, sep="/n")

# # make delimiter for heirarchical features | instead of ;
# # get vector of subjects from metadata that correspond to sample names in the feature table
# site_group_vec = metadata[match(rownames(community_df_t), metadata$sample_name), colnames(metadata)[METADATA_COL]]
# 
# # calculate indval index between taxa and each site group
# indval = multipatt(as.data.frame(community_df_t), site_group_vec)
# 
# # write out indicspecies results
# sink(paste(DATA_LOC, "otu_pipeline", "indicspecies_indicator_taxa.txt", sep="/"))
# cat(summary(indval))
# sink()
# 
# # approach 2) lefse
# # input data is created here, but lefse is run from the command line in otu_pipeline_script.txt
# 
# # append metadata grouping of interest to first row non-rarefied otu table for lefse
# lefse_community_df = rbind(site_group_vec, community_df)
# 
