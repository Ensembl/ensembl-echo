from ensembl.database import DBConnection
from ensembl.ncbi_taxonomy.api.utils import Taxonomy
from ensembl.ncbi_taxonomy.models import NCBITaxaNode, NCBITaxonomy
from collections import defaultdict
import csv 
import re
import sys

species = []
species_names = []
with open('/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/lepidoptera_species_unique_cores.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    #headings = next(reader)
    for row in reader:
        species.append((row[0].rstrip("\\n"),row[9]))
        species_names.append(row[0].rstrip("\\n"))

print("Number of input species", len(species))

unmatched = list(species)

ncbi_species = []
with open('/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/ncbi_species_subspecies_lepidoptera.csv', newline='') as cf:
    reader1 = csv.reader(cf,delimiter=',')
    for r in reader1:
        ncbi_species.append(r) 

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

matched_list = partial_matches(species,ncbi_species)
#print(len(matched_list))

def parse_cluster_file(filename):
    clusters = {}
    current_cluster = None

    with open(filename, "r") as f:
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

def between_left_right_index(session, parent_obj, tax_id):
    q = (session.query(NCBITaxaNode.taxon_id,NCBITaxaNode.parent_id, NCBITaxaNode.rank, NCBITaxaNode.left_index ,NCBITaxaNode.right_index, NCBITaxonomy.name)
        .join(NCBITaxonomy, NCBITaxaNode.taxon_id == NCBITaxonomy.taxon_id)
        .filter(NCBITaxaNode.left_index >= parent_obj.left_index, NCBITaxaNode.right_index <= parent_obj.right_index)
        .filter(NCBITaxonomy.name_class == 'scientific name')
        .filter(NCBITaxaNode.taxon_id != tax_id)
        .filter((NCBITaxaNode.rank == 'species') | (NCBITaxaNode.rank == 'subspecies')) #Species or subspecies
        .all()) 
    return q 

def get_parent_node(session, taxon_id):
    q = (session.query(NCBITaxaNode.parent_id)
            .filter(NCBITaxaNode.taxon_id == taxon_id))
    for i in q:
        #print(i[0])
        parent = Taxonomy.fetch_node_by_id(session, i[0])
        #print("parent",parent.taxon_id, "Left:", parent.left_index, "Right:", parent.right_index)

    return parent

def match_target_and_cluster_proteins(tax_id_of_protein, query_results, species_names, relative ):
    for row in query_results:
        if row[5] in species_names:
            #print("Tax id of proteins:", tax_id_of_protein)
            if  str(row[0]) in tax_id_of_protein:
                print(row)
                print("Tax id of the relatives:", row[0])
                if row[0] not in relative:
                    relative.append(row[0])


def get_closest_rel_within_cluster(number_of_relatives,matched_list, clusters, species_names, output_file):

    dbc = DBConnection('mysql://ensro@mysql-ens-meta-prod-1:4483/ncbi_taxonomy')
    result_dict = defaultdict(list)
    with dbc.session_scope() as session:
        for cluster_id, proteins in clusters.items():
            relative = []
            tax_id_of_protein = [p.split("_")[-1] for p in proteins]
            unique_tax_ids = list(dict.fromkeys(tax_id_of_protein))
            print("Cluster No.", cluster_id)
            print("Cluster members", proteins)
            print("Tax id of proteins:",tax_id_of_protein)
            if len(proteins) <= number_of_relatives or len(unique_tax_ids) <= number_of_relatives:
                print(proteins)
                #write code from get_fasta_of_closest_relatives.py to write these proteins to the fasta file
            else:       
                for tax_id, target in matched_list.items():
                    #print("Target genome:", target[0], target[2])
                    #taxon_id = target[0]
                    leaf = Taxonomy.is_leaf(session, tax_id)
                    #print(target[2])
                    print(tax_id)
                    print(target[0][1])
                    print("leaf",leaf)
                    if leaf == True:
                        parent = get_parent_node(session, tax_id)
                        desc = Taxonomy.num_descendants(session, parent.taxon_id)
                        print(desc)
                        query_results = between_left_right_index(session,parent, tax_id)
                        get_relatives = match_target_and_cluster_proteins(tax_id_of_protein, query_results, species_names, relative)
                        while len(relative) < number_of_relatives:
                            grandparent = get_parent_node(session, parent.taxon_id)
                            print("gp:", grandparent.taxon_id)
                            gp_query_results = between_left_right_index(session,grandparent, grandparent.taxon_id)
                            get_rel_gp = match_target_and_cluster_proteins(tax_id_of_protein, gp_query_results, species_names, relative)
                            parent = grandparent
                            #print("Before if condition", relative)
                            if Taxonomy.is_root(session, grandparent.taxon_id):
                                print("Root reached")
                                print("Relatives:", (cluster_id, relative))
                                result_dict[tuple(target)].append((cluster_id, tuple(relative)))
                                break
                        if len(relative) >= number_of_relatives:
                            print("Relatives:", (cluster_id, relative[:5]) )
                            result_dict[tuple(target)].append((cluster_id, tuple(relative[:5])))
                            break
                    else:
                        relative = []
                        desc = Taxonomy.num_descendants(session, tax_id)
                        print(desc)
                        parent = Taxonomy.fetch_node_by_id(session, tax_id)
                        query_result = between_left_right_index(session, parent, tax_id)
                        get_relative = match_target_and_cluster_proteins(tax_id_of_protein, query_result, species_names, relative)
                        while len(relative) < number_of_relatives:
                            grandparent = get_parent_node(session, parent.taxon_id)
                            print("gp:", grandparent.taxon_id)
                            gp_query_result = between_left_right_index(session,grandparent, grandparent.taxon_id)
                            get_rel_gp = match_target_and_cluster_proteins(tax_id_of_protein, gp_query_result, species_names, relative)
                            parent = grandparent
                            #print("Before if condition", relative)
                            if Taxonomy.is_root(session, grandparent.taxon_id):
                                print("Root reached")
                                print("Relatives:", (cluster_id, relative))
                                result_dict[tuple(target)].append((cluster_id, tuple(relative)))
                                break
                        if len(relative) >= number_of_relatives:
                            print("Relatives:", (cluster_id, relative[:5]) )
                            result_dict[tuple(target)].append((cluster_id, tuple(relative[:5])))
                            break
    print(result_dict)
    with open(output_file, 'w') as file:
        for target, data in result_dict.items():
            file.write(f"{target}: {data}\n")

    return result_dict


closest_relatives = get_closest_rel_within_cluster(5,matched_list, clusters, species_names, 'output_dict.txt' )
