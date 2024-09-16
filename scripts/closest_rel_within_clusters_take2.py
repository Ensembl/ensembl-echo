import pandas as pd 
from collections import defaultdict
import re
import requests
import time

def process_input_species_file(file_path):
    species = []
    #species_names = []
    df = pd.read_csv(file_path, delimiter='\t')
    for _, row in df.iterrows():
        species.append((row[0].rstrip("\\n"), row[9]))
        #species_names.append(row[0].rstrip("\\n"))
    
    return species

file_path = '/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/lepidoptera_species_unique_cores.csv'
species = process_input_species_file(file_path)

print("Number of input species", len(species))

def process_ncbi_species(ncbi_file_path):
    ncbi_species = []
    df = pd.read_csv(ncbi_file_path, delimiter=',')
    for _, row in df.iterrows():
        ncbi_species.append(row.tolist())

    return ncbi_species 

ncbi_file_path = '/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/ncbi_species_subspecies_lepidoptera.csv'
ncbi_species = process_ncbi_species(ncbi_file_path)

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
                matches.append((tax_id, sps_name,prod_name ))
                break

    for tax_id, sps_name, prod_name in matches:
        processed.setdefault(tax_id, []).append((sps_name, prod_name))

    return processed

matched_list = partial_matches(species, ncbi_species)
#print(matched_list)


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

clusters = parse_cluster_file("/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/results/rerun_with_unique_cores/lepidoptera_clusters_with_unique_cores.txt")
#print(clusters)

def get_rank(tax_id):

    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/dataset_report"
    response = requests.get(url)
    response.raise_for_status()  
    data = response.json()
    taxonomies = data.get('reports', [])
    for taxonomy in taxonomies:
        #rank_name = taxonomy.get('taxonomy', {}).get('classification', {}).get('rank', {}).get('name')
        rank = taxonomy.get('taxonomy', {}).get('rank')
    return  rank

def get_children(tax_id, max_retries=3, delay=5):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/dataset_report"
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
            data = response.json()
            taxonomies = data.get('reports', [])
            for taxonomy in taxonomies:
                children = taxonomy.get('taxonomy', {}).get('children')
                return children
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

def get_parents(tax_id, max_retries=3, delay=5):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/dataset_report"
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
            data = response.json()
            taxonomies = data.get('reports', [])
            for taxonomy in taxonomies:
                parents = taxonomy.get('taxonomy', {}).get('parents')
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


def related_ids(tax_id):
    url = url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/related_ids"
    response = requests.get(url)
    response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
    data = response.json()
    taxonomies = data.get('tax_ids', [])

    return taxonomies

def match_input_and_cluster_proteins(unique_tax_ids, desc, relative):
    for i in desc:
        #print(i)
        if i in unique_tax_ids and i not in relative:
            #print("Append to relative", i)
            relative.append(i)
            print("Append to relative", i)

def get_filtered_subtree(tax_id, max_retries=3, delay=5):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/filtered_subtree"
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
            data = response.json()
            taxonomies = data.get('reports', [])
            for taxonomy in taxonomies:
                children = taxonomy.get('taxonomy', {}).get('children')
                return children
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


def flatten(items):
    """Yield items from any nested iterable; see Reference."""
    for x in items:
        if isinstance(x, list) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x

def get_all_descendants(desc):
    all_desc = []
    for i in desc:
        tips = get_children(i)
        all_desc.append(tips)
    return list(flatten(all_desc))

def get_closest_rel_within_cluster(number_of_relatives,matched_list, clusters):

    #result_dict = defaultdict(list)
    for cluster_id, proteins in clusters.items():
        tax_id_of_protein = [int(p.split("_")[-1]) for p in proteins]
        unique_tax_ids = list(dict.fromkeys(tax_id_of_protein))
        print("Cluster No.", cluster_id)
        #print("Cluster members", proteins)
        #print("Tax id of the unique genomes the proteins belong to:",tax_id_of_protein)
        if len(unique_tax_ids) <= number_of_relatives:
            to_write = proteins
            #print(proteins)
            #write code from get_fasta_of_closest_relatives.py to write these proteins to the fasta file
        else:       
            for tax_id, target in matched_list.items():
                relative = []
                print("tax_id", tax_id, "target", target)
                #print(target[0][1])
                children = get_children(tax_id)
                print(f"Children of {tax_id}", children)
                parents = get_parents(tax_id)
                print(f"Parents of {tax_id}", parents)
                #if children == None:
                parent = parents[-1]
                #print(f"Immediate parent {parent}")
                desc = get_children(parent)
                print("descendants", desc)
                #print("Unique tax ids:", unique_tax_ids)
                get_relatives = match_input_and_cluster_proteins(unique_tax_ids, desc, relative)
                print("The tax ids of the genomes that are present in the cluster and in ensembl", relative)
                while len(relative) < number_of_relatives:
                    parents.pop()
                    gp = parents[-1]
                    print("gp", gp)
                    desc = get_children(gp)
                    print("descendants", desc)
                    tips = get_all_descendants(desc) #I need to get desc of tips,
                    #This logic is fine for the first gp, but when we go to the parent of this gp I will again miss the tips 
                    #because of which I am getting very few relatives. 
                    print("All descendents", tips)
                    get_gp_rel = match_input_and_cluster_proteins(unique_tax_ids, tips, relative)
                    print("The tax ids of the genomes that are present in the cluster and in ensembl (gp)", relative) 
                    if gp == 1:
                        print("Root reached")
                        get_rel_root = match_input_and_cluster_proteins(unique_tax_ids, tips)
                        #root_desc = get_children(gp)
                        #print("Closest relatives", (cluster_id, relative))
                        #print("The tax ids of the genomes that are present in the cluster and in ensembl (gp)", relative)
                        #print("Relatives:", target, ":", cluster_id, relative)
                        break
                if len(relative) >= number_of_relatives:
                    print("Relatives:", target, ":", cluster_id, relative)
                    break
                #print("relatives", relative)
    #print(result_dict)
    
    # with open(output_file, 'w') as file:
    #     for target, data in result_dict.items():
    #         file.write(f"{target}: {data}\n")

    return result_dict

#output_file = "results.csv"                           
closest_neighbours = get_closest_rel_within_cluster(5,matched_list, clusters)
#print(closest_neighbours)