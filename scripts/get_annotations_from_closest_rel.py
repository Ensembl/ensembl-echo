import os 
import csv
import subprocess

pep_files="/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/results/protein_seq_closest_sps_fasta_files"
core_with_tax_tsv="/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/core_names_tax_ids.tsv"

data = []
with open(core_with_tax_tsv, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    #headings = next(reader)
    for row in reader:
        data.append(row) 

target_genome_dir = "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/all_genome_files" 
results_dir = "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/results/all_gff_files"

for file in os.listdir(pep_files):
    filename = os.path.basename(file)
    tax_id = filename.split("_")[0]
    print(tax_id) 
    for core in data:
        if tax_id in core:
            #print(core)
            core_name = core[2]
            output_file = os.path.join(results_dir, f"{core_name}.gff")
            if not os.path.exists(output_file):
                target_genome = os.path.join(target_genome_dir, core[2] + ".fa")
                query_proteins = os.path.join(pep_files, file)
                command = f'{"miniprot"} -Iut16 --gff {target_genome} {query_proteins} > {output_file}'
                subprocess.run(command, shell=True)
            

