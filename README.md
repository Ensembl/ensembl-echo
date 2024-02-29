# Ensembl Conserved set of HOmologues (ECHO)

1. The script run_utils.py is run first to get the closest relatives of input core genomes. 
2. The script run_alfatclust.sh is then run to obtain clusters of input core peptide files.
3. Followed by the script get_fasta_of_closest_relatives.py, this script is used to get the fasta file for the target genome by selecting proteins of its closest relatives.
4. get_genome_files.sh is then run to obtain the genomes of the target species.
5. Finally, get_annotations_from_closest_relatives.py is used to run miniprot with the genomes obtained in step 4 to get gff files for the target genome based on proteins from its closest relatives. 
