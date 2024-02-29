#!/bin/bash


#Updated the dump_translation.pl file to include core name and taxonomy id in the fasta headers.
#Dumped the pep files using the dump_tanslations.pm script
# Combine the dumped files
cat *.fa >> 411_combined_lepidoptera_pep.fa

# Activate pyenv environment
pyenv activate ALFATClust

# Filter the pep file
filter_seqs -i 411_combined_lepidoptera_pep.fa -o 411_combined_lepidoptera_pep_filtered.fa

# Preprocess the pep file for alfatclust
replace_seq_header_spaces -i 411_combined_lepidoptera_pep_filtered.fa -o 411_combined_lepidoptera_pep_filtered_preprocessed.fa

# Run alfatclust with job submission
time ibsub -m 256gb alfatclust --evaluate lepidoptera_cluster_evaluation_201023.csv -i 411_combined_lepidoptera_pep_filtered_preprocessed.fa -o lepidoptera_clusters_with_updated_headers.txt
