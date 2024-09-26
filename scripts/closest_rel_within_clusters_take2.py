import pandas as pd
from collections import defaultdict
import re
import requests
import time
import pytaxonkit
import os


def process_input_species_file(file_path):
    species = []
    # species_names = []
    df = pd.read_csv(file_path, delimiter="\t")
    for _, row in df.iterrows():
        species.append((row[0].rstrip("\\n"), row[9]))
        # species_names.append(row[0].rstrip("\\n"))
    return species


def process_ncbi_species(ncbi_file_path):
    ncbi_species = []
    df = pd.read_csv(ncbi_file_path, delimiter=",")
    for _, row in df.iterrows():
        ncbi_species.append(row.tolist())
    return ncbi_species


def partial_matches(species, ncbi_species):
    matches = []
    processed = {}
    for sp in species:
        pattern = re.compile(re.escape(sp[0]))
        for sp2 in ncbi_species:
            if pattern.search(sp2[2]):
                tax_id = sp2[0]
                sps_name = sp2[2]
                prod_name = sp[1]
                matches.append((tax_id, sps_name, prod_name))
                break

    for tax_id, sps_name, prod_name in matches:
        processed.setdefault(tax_id, []).append((sps_name, prod_name))

    return processed


def parse_cluster_file(cluster_file):
    clusters = {}
    current_cluster = None

    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespaces
            if line.startswith("#Cluster"):
                # Extract cluster number
                cluster_number = line.split()[1]
                current_cluster = cluster_number
                clusters.setdefault(cluster_number, set())
            elif current_cluster and line:
                # Append line to current cluster
                clusters[current_cluster].add(line)

    # Convert sets to lists before returning
    for cluster_number, values_set in clusters.items():
        clusters[cluster_number] = list(values_set)

    return clusters

def get_parents_pytaxon(tax_id):
    result = pytaxonkit.lineage([tax_id])
    lineage_string = result["FullLineageTaxIDs"].tolist()
    for i in lineage_string:
        val = i.split(";")
        return list(map(int,val))

def get_parents(tax_id, max_retries=3, delay=5):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/dataset_report"
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
            data = response.json()
            taxonomies = data.get("reports", [])
            for taxonomy in taxonomies:
                parents = taxonomy.get("taxonomy", {}).get("parents")
                return parents
        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} of {max_retries} failed: {e}")
            if attempt < max_retries - 1:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
        except ValueError as e:
            print(f"An error occurred while processing the JSON response: {e}")
            break
        except KeyError as e:
            print(f"An error occurred while accessing data: {e}")
            break
    return None

def get_sequences_from_pep(protein_header, input_pep_files, outfile):
    
    split_header = protein_header.split("_")
    pep_file = '_'.join(split_header[2:8]) + "_translations.fa"
    prot_id = split_header[0]
    full_path = os.path.join(input_pep_files, pep_file)
    with open(full_path, "r") as file:
        write_sequence = False
        for line in file:
            if line.startswith(">"):
                if prot_id == line.split(" ")[0][1:]:
                    outfile.write(line)
                    write_sequence = True
                else:
                    write_sequence = False
            elif write_sequence:
                outfile.write(line)

def get_proteins_from_taxid(tax_id, proteins, target, output_path, input_pep_files):
    file_name = target + '_relatives.fa'
    file_path = os.path.join(output_path, file_name)
    if os.path.isfile(file_path):
        mode = "a"
    else:
        mode = "w"
    
    with open(file_path, mode) as file:
        for j in proteins:
            split_header = j.split("_")
            if str(tax_id) == split_header[-1]:
                print(j)
                seq = get_sequences_from_pep(j, input_pep_files, file)

def match_input_and_cluster_proteins(unique_tax_ids, desc, relative, number_of_relatives, proteins, target, output, input_pep_files):
    for i in desc:
        # print(i)
        if i in unique_tax_ids and i not in relative:
            relative.append(i)
            prot_ids = get_proteins_from_taxid(i, proteins, target[0][1], output, input_pep_files)
            if len(relative) == number_of_relatives:
                break

def get_all_desc(tax_id, desc):
    result = pytaxonkit.list([tax_id])
    for taxon, tree in result:
        subtaxa = [t for t in tree.traverse]
        for j in subtaxa:
            if j[1] == "species":
                desc.append(j[0])
    return desc

def get_lca(tax_ids):
    lca = pytaxonkit.lca(tax_ids)
    return lca

def get_rel_args(ncbi_file_path_arg, file_path_arg, clusters_arg):
    ncbi_file_path = ncbi_file_path_arg
    ncbi_species = process_ncbi_species(ncbi_file_path)
    file_path = file_path_arg
    species = process_input_species_file(file_path)
    print("Number of input species", len(species))
    matched_list = partial_matches(species, ncbi_species)
    clusters = parse_cluster_file(clusters_arg)
    lca = get_lca(list(matched_list.keys()))
    return matched_list, clusters, lca

def get_sequences_for_fasta(proteins, input_pep_files, output_dir):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "clusters_with_less_proteins.fa")
    
    with open(output_file, "a") as outfile:
        for i in proteins:
            seq = get_sequences_from_pep(i, input_pep_files, outfile)
                
def get_closest_rel_within_cluster(
    number_of_relatives, ncbi_file_path_arg, file_path_arg, clusters_arg, input_pep_files, output
):
    matched_list, clusters, lca = get_rel_args(
        ncbi_file_path_arg, file_path_arg, clusters_arg
    )
    for cluster_id, proteins in clusters.items():
        tax_id_of_protein = [int(p.split("_")[-1]) for p in proteins]
        unique_tax_ids = list(dict.fromkeys(tax_id_of_protein))
        print("Cluster No.", cluster_id)
        if len(unique_tax_ids) <= number_of_relatives:
            to_write = get_sequences_for_fasta(proteins, input_pep_files, output)
        else:
            for tax_id, target in matched_list.items():
                relative = []
                immi_desc = []
                print("tax_id", tax_id, "target", target)
                parents = get_parents_pytaxon(tax_id)
                print(f"Parents of {tax_id}", parents)
                print(type(parents))
                parents.pop()
                parent = parents[-1]
                immi_parent_desc = get_all_desc(parent, immi_desc)
                get_relatives = match_input_and_cluster_proteins(
                    unique_tax_ids, immi_desc, relative, number_of_relatives, proteins, target, output, input_pep_files
                )
                print(
                    "The tax ids of the genomes that are present in the cluster and in ensembl",
                    relative,
                )
                while len(relative) < number_of_relatives:
                    descendants = []
                    # Initialising the list in a loop because this will help reduce the search space.
                    # A great grand parent will have all the child nodes that a grand parent or a parent has.
                    parents.pop()
                    gp = parents[-1]
                    all_desc = get_all_desc(gp, descendants)
                    print(len(descendants))
                    get_gp_rel = match_input_and_cluster_proteins(
                        unique_tax_ids, descendants, relative, number_of_relatives, proteins, target, output, input_pep_files
                    )
                    print(
                        "The tax ids of the genomes that are present in the cluster and in ensembl (gp)",
                        relative,
                    )
                    if gp == lca:
                        print("LCA reached")
                        print(
                            "Relatives:",
                            target,
                            ":",
                            cluster_id,
                            relative,
                        )
                        
    return
