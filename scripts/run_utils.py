from ensembl.database import DBConnection
from ensembl.ncbi_taxonomy.api.utils import Taxonomy
from ensembl.ncbi_taxonomy.models import NCBITaxaNode, NCBITaxonomy
import csv
import re
import json

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
    for sp in species:
        pattern = re.compile(re.escape(sp[0]))
        for sp2 in ncbi_species:
            if pattern.search(sp2[2]):
                tax_id = sp2[0]
                sps_name = sp2[2]
                prod_name = sp[1]
                matches.append((tax_id, sps_name,prod_name ))
                break
    return matches

matched_list = partial_matches(species,ncbi_species)
#print(matched_list)
print("Number of species that match between input and ncbi api", len(matched_list))

def between_left_right_index(session, parent_obj, tax_id):
    q = (session.query(NCBITaxaNode.taxon_id,NCBITaxaNode.parent_id, NCBITaxaNode.rank, NCBITaxaNode.left_index ,NCBITaxaNode.right_index, NCBITaxonomy.name)
        .join(NCBITaxonomy, NCBITaxaNode.taxon_id == NCBITaxonomy.taxon_id)
        .filter(NCBITaxaNode.left_index >= parent_obj.left_index, NCBITaxaNode.right_index <= parent_obj.right_index)
        .filter(NCBITaxonomy.name_class == 'scientific name')
        .filter(NCBITaxaNode.taxon_id != tax_id)
        .filter((NCBITaxaNode.rank == 'species') | (NCBITaxaNode.rank == 'subspecies')) #Species or subspecies
        .all()) 
    return q 

def get_descendants(parent_obj ,desc, max_count, my_tax_id, session):
    if desc >= max_count:
        query_results = between_left_right_index(session,parent_obj, my_tax_id)
        num_desc_in_ensembl = 0
        #unique_num_desc_in_ensembl = 0
        unique_closest_relatives = []
        for row in query_results:
            if row[5] in species_names and row[0] not in unique_closest_relatives:
                print(row)
                unique_closest_relatives.append(row[0])
                num_desc_in_ensembl += 1
        print("Num desc in ensembl", num_desc_in_ensembl)
        #print(unique_closest_relatives)
        if len(unique_closest_relatives) > max_count:
            return unique_closest_relatives[:max_count]
        while num_desc_in_ensembl <  max_count:
            grandparent = get_parent_node(session, parent_obj.taxon_id)
            #grandparent = Taxonomy.parent(session, parent_obj.taxon_id)
            #print("grandparent",grandparent.taxon_id, "Left:", grandparent.left_index, "Right:", grandparent.right_index)
            query_results_gp = between_left_right_index(session,grandparent, my_tax_id)
            for row in query_results_gp:
                if row[5] in species_names and row[0] not in unique_closest_relatives:
                   print(row)
                   num_desc_in_ensembl += 1
                   unique_closest_relatives.append(row[0])
                   #unique_num_desc_in_ensembl += 1
            parent_obj = grandparent
            if Taxonomy.is_root(session, grandparent.taxon_id):
                break
        print("Num desc in ensembl", num_desc_in_ensembl)
        #print(unique_num_desc_in_ensembl)
        print(unique_closest_relatives)
        if num_desc_in_ensembl >= max_count:
            return unique_closest_relatives[:max_count]
    else:
        unique_closest_relatives = []
        num_desc_in_ensembl = 0
        #unique_num_desc_in_ensembl = 0
        while num_desc_in_ensembl < max_count:
            grandparent = get_parent_node(session, parent_obj.taxon_id)
            #grandparent = Taxonomy.parent(session, parent_obj.taxon_id)
            #print("grandparent",grandparent.taxon_id, "Left:", grandparent.left_index, "Right:", grandparent.right_index)
            query_results_gp = between_left_right_index(session,grandparent, my_tax_id)
            for row in query_results_gp:
                if row[5] in species_names and row[0] not in unique_closest_relatives:
                    print(row)
                    num_desc_in_ensembl += 1
                    unique_closest_relatives.append(row[0])
                    #unique_num_desc_in_ensembl += 1
            parent_obj = grandparent
            if Taxonomy.is_root(session, grandparent.taxon_id):
                break
        #print("Num of desc in ensembl", num_desc_in_ensembl)
        #print(unique_num_desc_in_ensembl)
        print(unique_closest_relatives)
        if len(unique_closest_relatives) > max_count:
            return unique_closest_relatives[:max_count]

        return unique_closest_relatives

def get_parent_node(session, taxon_id):
    q = (session.query(NCBITaxaNode.parent_id)
            .filter(NCBITaxaNode.taxon_id == taxon_id))
    for i in q:
        print(i[0])
        parent = Taxonomy.fetch_node_by_id(session, i[0])
        print("parent",parent.taxon_id, "Left:", parent.left_index, "Right:", parent.right_index)

    return parent

def get_closest_neighbours(matched_list, max_count, species_names, output_file_path):
    dbc = DBConnection('mysql://ensro@mysql-ens-meta-prod-1:4483/ncbi_taxonomy')
    with dbc.session_scope() as session:
        all_closest_rel_dict = {}
        for tup in matched_list:
            closest_rel_dict = dict()
            taxon_id = tup[0]
            leaf = Taxonomy.is_leaf(session, taxon_id)
            print(tup[2])
            print("leaf",leaf)
            if leaf == True:
                print("taxon_id",taxon_id)
                parent = get_parent_node(session, taxon_id)
                #parent = Taxonomy.parent(session, taxon_id)
                #print("parent",parent.taxon_id, "Left:", parent.left_index, "Right:", parent.right_index)
                desc = Taxonomy.num_descendants(session, parent.taxon_id)
                print(desc)
                #get_descendants(parent, desc, max_count,session)
                closest_rel_dict[taxon_id] = get_descendants(parent, desc, max_count, taxon_id, session)
            elif leaf == False:
                print("taxon_id",taxon_id)
                desc = Taxonomy.num_descendants(session, taxon_id)
                print(desc)
                node_by_id = Taxonomy.fetch_node_by_id(session, taxon_id)
                closest_rel_dict[taxon_id] = get_descendants(node_by_id, desc, max_count, taxon_id, session)
            #print(closest_rel_dict)
            all_closest_rel_dict.update(closest_rel_dict)
            print(all_closest_rel_dict)

    with open(output_file_path, 'w') as json_file:
        json.dump(all_closest_rel_dict, json_file)
            #json_file.write('\n')

    return closest_rel_dict
             
get_closest_neighbours(matched_list,10,species_names, "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/scripts/closest_relatives.json")
#dbc = DBConnection('mysql://ensro@mysql-ens-meta-prod-1:4483/ncbi_taxonomy')
#with dbc.session_scope() as session:
    #parent = Taxonomy.parent(session, "254363")
 #   p = get_parent_node(session, "254363")
    
  #  print(between_left_right_index(session, p, 254363))



