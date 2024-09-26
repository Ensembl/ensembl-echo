import argparse
from closest_rel_within_clusters_take2 import get_closest_rel_within_cluster

def main():
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('--num_of_rel', type=int, required=True)
    parser.add_argument('--ncbi_species_path', type=str, required=True)
    parser.add_argument('--species_path', type=str, required=True)
    parser.add_argument('--cluster_path', type=str, required=True)
    parser.add_argument('--input_pep_files', type=str, required=True)
    parser.add_argument('--output_dir', type=str, default="output")


    

    args = parser.parse_args()
    print(args)

    if(args.num_of_rel):
        output_file = "output.txt"                           
        closest_neighbours = get_closest_rel_within_cluster(args.num_of_rel, args.ncbi_species_path, args.species_path, args.cluster_path, args.input_pep_files, args.output_dir)
main()