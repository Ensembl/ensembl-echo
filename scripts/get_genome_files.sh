#!/bin/bash

csv_file="/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/lepidoptera_species_with_core.csv"

while IFS=',' read -r -a row; do
    column_data=${row[9]}
    target_file="${column_data}.fa"
    #echo "Checking $target_file..."
    #echo $column_data

    if [ ! -f "$target_file" ]; then
       echo "${column_data}"
       perl "${ENSEMBL_ROOT_DIR}/ensembl-analysis/scripts/sequence_dump.pl" \
            -host "mysql-ens-sta-5" \
            -port "4684" \
            -dbname $column_data \
            -user "ensro" \
            -nonref \
            --coord_system_name toplevel \
            -toplevel \
            -mask \
            -output_dir "/hps/nobackup/flicek/ensembl/compara/sbhurji/Development/ECHO_project/alfatclust_trial/data/all_genome_files/" \
            -filename "$target_file"
    else
         echo "BLAH"
    fi
done < $csv_file
