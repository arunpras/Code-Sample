
###########################################################################################
# Written by Arun Manoharan arun@primediscoveries.com 6/12/18
# This script analyses an otu table and produces data for visualizations
###########################################################################################

### TO-DO

# COMMAND LINE ARGUMENTS
# take in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# test if correct number of arguments, if not return an error
NUM_ARGS <- 2
if (length(args) != NUM_ARGS) {
  stop(paste("fastq_qc_visualization.R requires", NUM_ARGS, "arguments:
       PROJECT_PATH: absolute path to project directory
       DATABASE: database that taxonomy assignment was based on; options are ['gg_13_8', 'hitdb', 'silva']

       This script requires files produced by otu_pipeline_script.txt in PROJECT_PATH/otu_pipeline
       
       Example Usage:
       Rscript otu_visualization.R
       /home/arun/workdir/test_datasets/qiita_single_fastq_test
       'gg_13_8'", sep=" "), call. = F)
}

# REQUIREMENTS

# requirements for phyloseq object creation
### currently, phyloseq is only installed in R in qiime2-2018.4 environment on ec2
library(phyloseq)
library(readr)
library(reshape2)

# EXAMPLE USAGE

### currently, phyloseq is only installed in R in qiime2-2018.4 environment on ec2
# conda activate /home/ubuntu/miniconda2/envs/qiime2-2018.4

# Rscript otu_visualization.R \
# /home/arun/workdir/test_datasets/qiita_single_fastq_test \
# "gg_13_8" \

###########################################################################################
# SETUP - ARGUMENTS, DIRECTORY CREATION

 

# assign variables based on command line arguments
PROJECT_PATH <- args[1] # absolute path to project directory
DATABASE <- args[2] # database that taxonomy assignment was based on; options are ['gg_13_8', 'hitdb', 'silva']

#PROJECT_PATH <- '/Volumes/prime/microbiome_downloads/psoriasis/NCBI/PRJEB25915' # absolute path to project directory
#DATABASE <- 'gg_13_8' # database that taxonomy assignment was based on; options are ['gg_13_8', 'hitdb', 'silva']

# make directories for output
dir.create(paste(PROJECT_PATH, "otu_visualization_data", sep="/"))
OUTPUT_LOC <- paste(PROJECT_PATH, "otu_visualization_data", sep="/")
dir.create(paste(OUTPUT_LOC, "1_tree", sep="/"))
dir.create(paste(OUTPUT_LOC, "2_abundance", sep="/"))
dir.create(paste(OUTPUT_LOC, "3_b_f_ratio", sep="/"))
dir.create(paste(OUTPUT_LOC, "other_files", sep="/"))
cat(paste(Sys.time(), "Creating directories for output completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"))

# set paths for data files produced by otu_pipeline_script.sh
FEATURE_TABLE_PATH <- paste(PROJECT_PATH, "otu_pipeline", "1_denoise", "feature_table", sep="/")
TAXONOMY_PATH <- paste(PROJECT_PATH, "otu_pipeline", "3_taxonomy", sep="/")
PHYLOGENY_PATH <- paste(PROJECT_PATH, "otu_pipeline", "2_phylogeny", sep="/")

###########################################################################################
# GET TREE FILE IN CORRECT DIRECTORY
cat(paste(Sys.time(), "Moving tree and taxonomy files to otu_data_visualization directory started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
print("copying tree file from otu_pipeline directory to otu_data_visualization 1_tree directory")
system(paste("cp ", PHYLOGENY_PATH, "/tree.nwk ", OUTPUT_LOC, "/1_tree/asv_tree.nwk", sep=""))

# GET TAXONOMY MAPPING FILE IN CORRECT DIRECTORY
print("copying taxonomy-ASV code mapping file from otu_pipeline directory to otu_data_visualization directory")
system(paste("cp ", TAXONOMY_PATH, "/taxonomy.tsv ", OUTPUT_LOC, "/other_files/taxonomy.tsv", sep=""))
cat(paste(Sys.time(), "Moving tree and taxonomy files to otu_data_visualization directory completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")

###########################################################################################
# CREATE PHYLOSEQ OBJECT
cat(paste(Sys.time(), "Creating phyloseq object started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
print("creating phyloseq object")


# load data - feature table
cat(paste(Sys.time(), "Loading feature table started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
qiime_data = import_biom(BIOMfilename = paste(FEATURE_TABLE_PATH, "feature-table_json.biom", sep="/"))
# get otu table as data frame
OTU = otu_table(qiime_data)
cat(paste(Sys.time(), "Loading feature table completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")


# load data - taxonomy table
cat(paste(Sys.time(), "Loading taxonomy table started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
taxonomy = read.csv(paste(TAXONOMY_PATH, "taxonomy.tsv", sep="/"), sep = "\t")
# for greengenes, taxonomy follows standard format k, p, c, o, f, g, s as in gg_13_8 99% otus the taxonomy classifier was trained on
if (DATABASE == 'gg_13_8') {
  taxmat <- as.matrix(colsplit(taxonomy$Taxon, pattern = "; ", names = c("k", "p", "c", "o", "f", "g", "s")))

# for hitdb, taxonomy follows standard format p, c, o, f, g, s as in HITdb_taxonomy_qiime.txt the taxonomy classifier was trained on
# kingdom column must be manually added for hitdb (the only archaea belong to the phylum Euryarchaeota)
} else if (DATABASE == 'hitdb') {
  tax_df <- colsplit(taxonomy$Taxon, pattern = ";", names = c("p", "c", "o", "f", "g", "s"))
  # assign kingdom: for every line in taxmat, determine if kingdom is Bacteria (phylum != Euryarchaeota) or Archaea (phylum == Euryarchaeota)
  tax_df$k = ""
  for (i in seq(1, dim(tax_df)[1], 1)) {
    if (tax_df[i, "p"] == "Unassigned") {
      tax_df[i, "k"] <- "Unassigned"
      tax_df[i, "p"] <- ""
    } else if (tax_df[i, "p"] == 'Euryarchaeota') {
      tax_df[i, "k"] <- "Archaea"
    } else {
      tax_df[i, "k"] <- 'Bacteria'
    }
  }
taxmat <- as.matrix(tax_df)

# for silva, tax table must be coerced to greengenes format using crossclassify
} else if (DATABASE == 'silva') {
  print("translation from silva to greengenes taxonomy standard naming is not implemented yet")
}

# format taxonomy matrix with OTU feature IDs as rownames, standardized taxonomic levels as column names
rownames(taxmat) <- taxonomy$Feature.ID
colnames(taxmat) <- c("k", "p", "c", "o", "f", "g", "s")

# create phyloseq taxonomy table object from taxonomy matrix
TAX = tax_table(taxmat)
cat(paste(Sys.time(), "Loading taxonomy table completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")

# load data - tree
TREE <- read_tree(paste(OUTPUT_LOC, "1_tree", "asv_tree.nwk", sep="/"))

# create phloseq object
ps_object = merge_phyloseq(OTU, TAX, TREE)
# save phyloseq object
saveRDS(ps_object, file = paste(OUTPUT_LOC, "phyloseq_object.RDS", sep="/"))

cat(paste(Sys.time(), "Creating phyloseq object completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")

###########################################################################################
# EXPORT TREE WITH TAXA AS TIP LABELS
cat(paste(Sys.time(), "Exporting phylogentic tree with ASV codes and taxa names as tip names started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
print("exporting phylogenetic tree")

# export tree with taxa names (full lineage)
tree <- phy_tree(ps_object)
tree$tip.label <- taxonomy$Taxon[match(tree$tip.label, taxonomy$Feature.ID)]
ape::write.tree(tree, paste(OUTPUT_LOC, "1_tree", "taxa_tree.nwk", sep="/"))

# export tree with only top TOP_X_PERCENT % of ASVs
TOP_X <- 50
# prune less abundant ASVs from phyloseq object
top_asvs = names(sort(taxa_sums(ps_object), decreasing = TRUE)[1:min(TOP_X, dim(OTU)[1])])
ps_object_pruned = prune_taxa(top_asvs, ps_object)
# export tree
tree_pruned <- phy_tree(ps_object_pruned)
tree_pruned$tip.label <- taxonomy$Taxon[match(tree_pruned$tip.label, taxonomy$Feature.ID)]
ape::write.tree(tree_pruned, paste(OUTPUT_LOC, "1_tree", paste("taxa_tree_pruned_to_", TOP_X,".nwk", sep=""), sep="/"))

# export tree with tip labels as either "genus species" or "<finest taxonomic rank> other"
tree_pruned <- phy_tree(ps_object_pruned)

# generate labels for tree from Taxon column of taxonomy file split on taxa delimiter (tax_df)
taxonomy_split <- data.frame(taxonomy, colsplit(taxonomy$Taxon, pattern = ";", names = c("k", "p", "c", "o", "f", "g", "s")))

# initialzize column for tree label
taxonomy_split$tree_label <- ""

taxa_levels <- c("k", "p", "c", "o", "f", "g", "s")
taxa_level_names <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# for every row of taxonomy_split file, make label "genus species" if those are classified, otherwise make label "<rank> <name>" at finest level possible
for (i in seq(1, dim(taxonomy_split)[1])) {
  # get species column entry from taxonomy_split from row i
  species <- taxonomy_split[i, "s"]
  # if species classification exists, make tree label "genus species"
  species <- strsplit(species, "s__")[[1]][2]
  if (!is.na(species)) {
    # eliminate leading "s__" tag from species name
    genus <- strsplit(taxonomy_split[i, "g"], "g__")[[1]][2]
    taxonomy_split[i, "tree_label"] <- paste(genus, species, sep=" ")
  # if species label does not exist, make tree label based on finest taxonomic level that is classified
  } else {
    j <- length(taxa_levels) - 1
    # while taxa are undefined, continue moving to coarser taxonomic level
    while ((is.na(taxonomy_split[i, taxa_levels[j]]) || taxonomy_split[i, taxa_levels[j]]=="")&& j>1) {
      j <- j - 1
    }
    # make label "<rank> <name>" at finest level possible
    if(j==0){
      print(taxonomy_split[i, "tree_label"] <- paste(taxa_level_names[1], strsplit(taxonomy_split[i, taxa_levels[1]], paste(taxa_levels[1], "__", sep=""))[[1]][1], sep=" "))
      taxonomy_split[i, "tree_label"] <- paste(taxa_level_names[1], strsplit(taxonomy_split[i, taxa_levels[1]], paste(taxa_levels[1], "__", sep=""))[[1]][1], sep=" ")
    }else{
    # make label "<rank> <name>" at finest level possible
    taxonomy_split[i, "tree_label"] <- paste(taxa_level_names[j], strsplit(taxonomy_split[i, taxa_levels[j]], paste(taxa_levels[j], "__", sep=""))[[1]][2], sep=" ")
    }
  }
}

tree_pruned$tip.label <- taxonomy_split$tree_label[match(tree_pruned$tip.label, taxonomy_split$Feature.ID)]
ape::write.tree(tree_pruned, paste(OUTPUT_LOC, "1_tree", paste("taxa_tree_pruned_to_", TOP_X,".nwk", sep=""), sep="/"))
cat(paste(Sys.time(), "Exporting phylogentic tree with ASV codes and taxa names as tip names completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")

###########################################################################################
# GET TAXONOMIC DISTRIBUTION
cat(paste(Sys.time(), "Generating abundance tables started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
print("analysing abundance")

# write otu abundance data in heirarchical, uncollapsed form
#names(df)<-make.names(namesdf))

# get otu table (indexed by Feature ID)
otu_tab <- otu_table(ps_object)
# get taxonomy assignment table (indexed by Feature ID)
tax_tab <- tax_table(ps_object)
# merge otu and taxonomy information based on Feature IDs
heirarchical_abundance <- merge(otu_tab, tax_tab, by.x = 0, by.y = 0)
#heirarchical_abundance <- cbind(otu_tab, tax_tab) ##changed this to fix error thrown
# change Row.names column name to Feature_ID
colnames(heirarchical_abundance)[[1]] <- "Feature_ID"
# write heirarchical data out to file
write.table(heirarchical_abundance, file = paste(OUTPUT_LOC, "2_abundance", paste("absolute_abundance", ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)

if (DATABASE != 'silva') {

  # summarize otu table by taxonomic level of interest
  # NArm = T ignores any ASV for which taxonomy assignment is zero, so relative abundance of all categorized taxa sums to 1
  data_tax_glom_k <- tax_glom(ps_object, 'k', NArm = F)
  data_tax_glom_p <- tax_glom(ps_object, 'p', NArm = F)
  data_tax_glom_c <- tax_glom(ps_object, 'c', NArm = F)
  data_tax_glom_o <- tax_glom(ps_object, 'o', NArm = F)
  data_tax_glom_f <- tax_glom(ps_object, 'f', NArm = F)
  data_tax_glom_g <- tax_glom(ps_object, 'g', NArm = F)
  data_tax_glom_s <- tax_glom(ps_object, 's', NArm = F)

  # transform otu table to relative abundance
  data_tax_glom_rarefied_k <- transform_sample_counts(data_tax_glom_k, function(OTU) OTU/sum(OTU) )
  data_tax_glom_rarefied_p <- transform_sample_counts(data_tax_glom_p, function(OTU) OTU/sum(OTU) )
  data_tax_glom_rarefied_c <- transform_sample_counts(data_tax_glom_c, function(OTU) OTU/sum(OTU) )
  data_tax_glom_rarefied_o <- transform_sample_counts(data_tax_glom_o, function(OTU) OTU/sum(OTU) )
  data_tax_glom_rarefied_f <- transform_sample_counts(data_tax_glom_f, function(OTU) OTU/sum(OTU) )
  data_tax_glom_rarefied_g <- transform_sample_counts(data_tax_glom_g, function(OTU) OTU/sum(OTU) )
  data_tax_glom_rarefied_s <- transform_sample_counts(data_tax_glom_s, function(OTU) OTU/sum(OTU) )

  # join taxonomic information to otu table
  tax_plot_data_k = merge(otu_table(data_tax_glom_rarefied_k), tax_table(data_tax_glom_rarefied_k)[, 'k'], by.x=0, by.y=0)
  tax_plot_data_p = merge(otu_table(data_tax_glom_rarefied_p), tax_table(data_tax_glom_rarefied_p)[, 'p'], by.x=0, by.y=0)
  tax_plot_data_c = merge(otu_table(data_tax_glom_rarefied_c), tax_table(data_tax_glom_rarefied_c)[, 'c'], by.x=0, by.y=0)
  tax_plot_data_o = merge(otu_table(data_tax_glom_rarefied_o), tax_table(data_tax_glom_rarefied_o)[, 'o'], by.x=0, by.y=0)
  tax_plot_data_f = merge(otu_table(data_tax_glom_rarefied_f), tax_table(data_tax_glom_rarefied_f)[, 'f'], by.x=0, by.y=0)
  tax_plot_data_g = merge(otu_table(data_tax_glom_rarefied_g), tax_table(data_tax_glom_rarefied_g)[, 'g'], by.x=0, by.y=0)
  tax_plot_data_s = merge(otu_table(data_tax_glom_rarefied_s), tax_table(data_tax_glom_rarefied_s)[, 's'], by.x=0, by.y=0)

  # rename "Row.names" column to Feature ID
  colnames(tax_plot_data_k)[[1]] <- "Feature_ID"
  colnames(tax_plot_data_p)[[1]] <- "Feature_ID"
  colnames(tax_plot_data_c)[[1]] <- "Feature_ID"
  colnames(tax_plot_data_o)[[1]] <- "Feature_ID"
  colnames(tax_plot_data_f)[[1]] <- "Feature_ID"
  colnames(tax_plot_data_g)[[1]] <- "Feature_ID"
  colnames(tax_plot_data_s)[[1]] <- "Feature_ID"

  # write out taxonomic information for each level
  write.table(tax_plot_data_k, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 'k', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  write.table(tax_plot_data_p, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 'p', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  write.table(tax_plot_data_c, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 'c', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  write.table(tax_plot_data_o, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 'o', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  write.table(tax_plot_data_f, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 'f', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  write.table(tax_plot_data_g, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 'g', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  write.table(tax_plot_data_s, file = paste(OUTPUT_LOC, "2_abundance", paste("percent_abundance_", 's', ".tsv", sep=""), sep="/"), sep='\t', quote = F, row.names = F)
  cat(paste(Sys.time(), "Generating abundance tables completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
  
  # write out Bacteroidetes:Firmicutes ratio
  cat(paste(Sys.time(), "Writing out Bacteroidetes:Firmicutes ratio started.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")

  # get rows of OTU table (with taxonomic information appended, aggregated at the phylum level) corresponding to Bacteroidetes and Firmicutes
  brow <- tax_plot_data_p[which(tax_plot_data_p$p == "p__Bacteroidetes"), ]
  frow <- tax_plot_data_p[which(tax_plot_data_p$p == "p__Firmicutes"), ]
  # limit rows to samples only, exclude OTU code (Row.names) and phylum (p) columns
  brow <- subset(brow, select = -c(Feature_ID, p))
  frow <- subset(frow, select = -c(Feature_ID, p))
  # divide relative abundance of Bacteroidetes in each sample by relative abundance of Firmicutes and transpose from row vector to column vector
  bf_ratio = t(brow/frow)
  # rename column vector header to B_F_ratio
  colnames(bf_ratio) <- "3_b_f_ratio"
  # add sample name column
  bf_ratio <- data.frame("sample_name"=rownames(bf_ratio), bf_ratio)
  # write out B:F ratio to file
  write.table(bf_ratio, file = paste(OUTPUT_LOC, "3_b_f_ratio", "b_f_ratio.tsv", sep="/"), sep='\t', quote=F, row.names=F)
  cat(paste(Sys.time(), "Writing out Bacteroidetes:Firmicutes ratio completed.\n", sep=" "), file=paste(OUTPUT_LOC, "otu_visualization_log.txt", sep="/"), append=T, sep="/n")
} else {
  print("silva taxonomy mapping to k, p, c, o, f, g, s not implemented yet")
}
