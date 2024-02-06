# Steps for generating metadata:

## select_samples.py
Aligns kmers with its corresponding metadata. Specify specific the metadata column you want as well as the specific labels. If you want all possible labels, pass `all` to the `--desired_features` argument

### Parameters:
-`--kmer_file`:filepath to kmers
-`--kmer_meta`: filepath to kmer metadata
-`--sample_meta`: filepath to sample metadata
-`--meta_col`: metadata column you want
-`--desired_features`: features you want seperated by columns, use 'all' for all features
-`--meta_save_dest`: destination to save metadata to
-`--kmer_save_dest`: destination to save kmers to

### example:
The following script selects kmers with metadata label in 'UBERON:feces,UBERON:skin,UBERON:nostril'
Command: 
```bash
python /Users/prime_galen/Documents/bioinformatics/scripts/data_curation/select_samples.py \
	--kmer_file '/Users/prime_galen/Documents/bioinformatics/test_datasets/449/449_6-mers_1000.npz' \
	--kmer_meta '/Users/prime_galen/Documents/bioinformatics/test_datasets/449/449_6-mers_1000_meta'\
	--sample_meta '/Users/prime_galen/Documents/bioinformatics/test_datasets/449/449_metadata.txt' \
	--meta_save_dest '/Users/prime_galen/Documents/bioinformatics/test_datasets/449/test_metadata.csv'\
	--kmer_save_dest '/Users/prime_galen/Documents/bioinformatics/test_datasets/449/test_metadata.npz'\
	--meta_col 'body_habitat' \
	--desired_features 'UBERON:feces,UBERON:skin,UBERON:nostril' \
```
Output: 
```bash
Selecting following labels from: body_habitat
UBERON:feces
UBERON:skin
UBERON:nostril
Saved 448 kmers
```





