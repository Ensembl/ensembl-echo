import csv
import re
import os
import json

cluster_csv = "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/results/rerun_with_updated_headers/lepidoptera_cluster_evaluation_231023.csv"

centre_seqs = []
with open(cluster_csv, newline='') as clust_eval:
     reader2 = csv.reader(clust_eval, delimiter=',')
     headings = next(reader2)
     for clusters in reader2:
         centre_seq = clusters[4]
         centre_seq_tax_id = centre_seq.split("_")[-1]
         centre_seqs.append(clusters[4])


#closest_relative_dict = {'753202': [214277, 987995, 320037, 987925, 988049, 997545, 1857961, 214171, 875885, 987895], '1594315': [989769, 1035111, 1100915, 1100916, 1100963, 1101027, 1869985, 753214, 758706, 758717], '55057': [271217, 875884, 987933, 987983, 987985, 988041, 988125, 997540, 116126, 116130], '938226': [987859, 55057, 116126, 116130, 179674, 214171, 214277, 254363, 271217, 320037], '721163': [721137, 721165, 7116, 7130, 33412, 33443, 33448, 42275, 55057, 64459]}

def get_proteins_of_closest_species(json_file_path, centre_seqs, dir_path_pep_files, output_dir):
    with open(json_file_path, 'r') as json_file:
        closest_relative_dict = json.load(json_file)
    
    for k,v in closest_relative_dict.items():
        print(k, ":", v)
        all_pep_seqs = []
        for cs in centre_seqs:
            centre_seq_tax_id = cs.split("_")[-1]
            if centre_seq_tax_id in str(v):
                core_name = cs.split("_")[2:-1]
                pep_file_name = "_".join(core_name) + "_translations.fa"
                #print(pep_file_name)
                full_path = os.path.join(dir_path_pep_files, pep_file_name)
                with open(full_path, "r") as file:
                    printing = False
                    pep_seqs = []
                    for line in file:
                        if line.startswith(">"):
                            protein_id_cs = cs.split("_")[0]
                            protein_id_pep_file = line.split()[0].strip(">")
                            printing = False
                            if protein_id_cs == protein_id_pep_file:
                                printing = True 
                                pep_seqs.append(line.rstrip() + '\n')
                                #print(k)
                                #print(cs)
                                #print(line.rstrip())
                        elif printing:
                            #print(line.rstrip())
                            pep_seqs[-1] += line.rstrip()
                            #pep_seqs.append(line.rstrip())
                    all_pep_seqs.extend(pep_seqs)

        output_file_path = os.path.join(output_dir, f'{k}_selected_proteins.fasta')
        with open(output_file_path, 'w') as output_file:
            output_file.write('\n'.join(all_pep_seqs) + '\n')
        print(f"File written: {output_file_path}")

json_file_path = "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/scripts/closest_relatives.json"
dir_path_pep_files = "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/all_pep_files"
output_dir = "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/results/protein_seq_closest_sps_fasta_files/"
get_proteins_of_closest_species(json_file_path, centre_seqs,dir_path_pep_files, output_dir)
