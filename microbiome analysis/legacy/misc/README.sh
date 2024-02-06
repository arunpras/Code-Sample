############################################################################################################################
'''
1. Script to process fastq
creates a q2_importable folder with sequences
'''
cd <fastq folder>
source activate qiime1
bash /home/ubuntu/microbiome_downloads/make_importable.sh *.fastq

#remove part of filename from the files
cd q2_importable
rename 's/.fastq_/_/g' *.fastq.gz
sed -i -e 's/.fastq_/_/g' MANIFEST

#this creates usable folders with sequences
############################################################################################################################
'''
2. Generating Qiime outputs

make_biom.sh:
takes the fastq.gz in q2_importable, preprocesses, and trims at 100

location: 
	/home/arun/bioinformatics/arun/preprocessing/make_biom.sh

demo script: /home/arun/bioinformatics/arun/1070/commands.sh

'''
usage: 
 	./make_biom <q2_importable path> <metadata path> <dataset name> <output folder>

example usage:
	cd /home/arun/bioinformatics/arun/preprocessing
	./make_biom.sh ../1070/q2_importable ../1070/1070_metadata.txt 1070 ../1070	

/home/arun/bioinformatics/preprocessing/generate_OTU.sh we 449_metadata.txt 449 blha

output:
	dataset_name_input/
	    dataset_name_demux.qzv                      qimme demux summarize output  
	    dataset_name_input_data.qza                 qiime object after intial data import  
	biom/  
	    feature-table.biom                          OTU feature table saved as a BIOM file  
	deblur/  
	    post_deblur/  
	            deblur-stats.qza                    qiime deblur denoise-16S output: per sample stats  
	            post_deblur_rep_seq_summary.qzv     qiime feature-table tabulate-seqs output: provides mapping of each feature to its sequence
	            post_deblur_summarize_table.qzv     qiime feature-table summarize output: provides info on how many sequences associated with each feature
	            representative_sequences.qza        qiime deblur denoise-16S output: resulting feature sequences  
	            table.qza                           qiime deblur denoise-16S output: resulting denoised feature table  
	    qscores/  
	            filtered_sequences.qza              qiime quality-filter q-score output, filtered sequences based on quality scores  
	            filter_stats.qza                    qiime quality-filter q-score output, stats from quality filter  
	deblur.log                                      deblur log  
	post_qiime_processing/                           
	    dataset_name_FINAL_deblur_feature.txt       OTU frequency feature table, rows are samples, columns as OTUs  

############################################################################################################################
'''
3. Generating OTU feature and label tables:

merge_features_metadata.R:
Takes in OTU feature tsv from qiime2 and a metadata.txt. Aligns samples between feature table and metadata and saves as X_features and Y_metadata

location: 
	/home/arun/bioinformatics/arun/preprocessing/merge_features_metadata.R
'''
usage:
	Rscript <feature table> <metadata.txt> <output name>

example usage:
	cd /home/arun/bioinformatics/arun/preprocessing
	Rscript ../preprocessing/merge_features_metadata.R post_qiime_processing/1070_FINAL_deblur_feature.txt 1070_metadata.txt classification/

output:
	classification/X_features.csv 	#OTU feature table/Rows are samples, Columns is OTU freq
	classification/Y_metadata.txt 	#Metadata/Rows are samples, Columns are Metadata
############################################################################################################################
'''
4. Kmer generation:
micropheno.py --genkmer

location:
	/home/arun/bioinformatics/deep_learning/Micropheno
'''

#check if q2 importable needs to be unzipped
usage:
	python3 micropheno.py 
		--genkmer 		#generate kmer
		--inaddr 		#path to files
		--out 			#output folder
		--filetype 		#fasta/fastq/fsa 
		--cores 		#number of cores 
		--KN 			#number of k:subsampling depth ex: 6:1000 (6-mer with 1000 samples)
		--name 			#name of output


###FILTER KMER
###WHICH ones to remove
#FASTQC trimmer
#OPTIMIZE Sampling depth/bootstrap
example usage:
	python3 /home/arun/bioinformatics/deep_learning/Micropheno/micropheno.py --genkmer --inaddr /home/ubuntu/microbiome_downloads/449/q2_importable --out 449_kmer/kmers/ --filetype fastq --cores 10 --KN 6:1000 --name 449


python3 micropheno.py --genkmer --inaddr /Users/Prime/Downloads/ssr_83462 --out /Users/Prime/Downloads/ssr_83462/kmers --filetype fastq --cores 10 --KN 6:1000 --name eugene

output:
	449_6-mers_1000_log 	#log file
	449_6-mers_1000_meta	#metadata
	449_6-mers_1000.npz		#actual kmers (basically just a pickled sparse matrix)

#add commands to open and edit files (npz)


############################################################################################################################
'''
5. Running kmer through MicroPheno neural net:
micropheno.py --train_predictor

location:
	/home/arun/MicroPheno
'''
Summary of steps:
1. Generate kmers
2. Generate Y file
3. Run through DNN
4. Visualize results
5. Ordination
###
'''
Combine Kmers

location:
	/home/arun/MicroPheno/combine_kmers.py
'''
usage:
	python combine_kmers,py
		--kmer_1			#first kmer file
		--kmer_2			#additional kmers to add to kmer_1 (kmers must both need to be same sampling depth and the same)
		--num_add_x			#how many samples from kmer_1 to add to kmer_2. Useful when you only want a select number from an additional dataset
		--filename			#name of output file
	

example usage:
	cd /home/arun/MicroPheno/test_kmer_combine
	python ../combine_kmers.py --kmer_1='../449_kmer/kmers/449_6-mers_1000.npz' --kmer_2='../1939_kmer/1939_6-mers_1000.npz' --num_add_x='200' --filename='449_1939_200_6_1000'

output:
	449_1939_200_6_1000.npz	#combination of 449 and 1939 datasets, with only the top 200 from 1939, 6mer, sampling depth of 1000
###
'''
Training the neural net

location:
'''
usage:
	python3 micropheno.py 
		--train_predictor 	#train a predictor
		--model DNN 		#use a deep neural net
		--arch 				#DNN architecture in following format: num_neurons,dropout rate (btwn 0 and 1)
							#ex: 1024,0.2,512,0.1,128,64 
		--batchsize 		#training batch size
		--epochs  			#num_epochs
		--x 				#kmer file
		--y 				#single column text file, with label per line (see above to generate y file)
		--name 				#name 
		--out 				#output directory

example usage:
	python3 /home/arun/bioinformatics/deep_learning/Micropheno/micropheno.py --train_predictor --model DNN --arch 1024,0.2,512,0.1,128,64 --batchsize 10 --epochs 100 --x 449_6-mers_1000.npz --y 449_bodysite_Y_2removed.txt  --name 449  --out output_dir/

output:
	output_dir/nn_classification_results_449_layers_mlp_2048-0.2-1024-0.2-512-0.1-128-64_0.59.pickle 	#pickled file of the layers
	output_dir/nn_classification_results_449_mlp_2048-0.2-1024-0.2-512-0.1-128-64_0.59.pickle			#pickled file of the results
############################################################################################################################

#make sure all the pipelines work
(take it from qiita fastq, single file with q2 importable)
another from the ncbi
(download one from qiita and one from ncbi, script should run on both) (generate out taxonomy)

Every folder should have a readme
