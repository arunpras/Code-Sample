## Written By Arun Manoharan
## Kmer Scripts
### bootstrapping
Example script:
```
python3 micropheno.py --bootstrap --indir /home/galen/bioinformatics/to_be_moved/449/q2_importable_nongz --out output_dir/ --filetype fastq --kvals 3,4,5,6 --nvals 10,100,200 --name 449
```
### generate_kmers.py
This script generates Kmers

Arguments:\
--inaddr 		input_folder\
--out 			output folder\
--filetype 		(fasta/fastq/fsa/fastq.gz)\
--cores 		num cores to use\
--KN 			k:sampling depth\
--name 			name of output

Example usage:
```
python generate_kmers.py --inaddr ../../../projects/Bharath/q2 --out kmers/ --filetype fastq.gz --cores 10 --KN 4:2000,6:5000 --name arun
```

### train_DNN.py
Trains and tests a NN from the givin input
Arguments:\
--model DNN\
--arch model architecture\
--batchsize batch_size\
--epochs num_epochs\
--x kmer_file\
--y metadata\
--name name of output\
--out output_director

Example usage:
```
python /Users/arun/Documents/bioinformatics/scripts/kmer_processing/Micropheno/train_DNN.py --model DNN --arch 1024,0.2,512,0.1,128,64 --batchsize 10 --epochs 50 --x kmers/449_b_e_6_1000.npz --y 449_b_e_bodyhabitat.txt  --name 449_b_e  --out output_dir/
```

### visualization.py
Runs a TSNE and outputs the ordination

All the inputs are in the yaml file, see projects/449_b_e_bodyhabitat/440_bharath_eugene.yaml

Example usage:
```
python /Users/arun/Documents/bioinformatics/scripts/kmer_processing/Micropheno/visualization.py 449_bharath_eugene.yaml
```