#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// db connection
params.db = ''
params.host = 'mysql-ens-sta-5'
params.port = '4684'
params.user = 'ensro'
params.pass = ''
params.csvFile = ''
params.dump_params = '--canonical_only'

//repos
params.enscode = ''

csvFile = file(params.csvFile)

/*
process PER_GENOME {
    input:
        val genome
    output:
        path "out_path/*.txt", emit: txt
    script:
    """
    mkdir -p out_path
    echo ${genome} > out_path/${genome}.txt;
    """
}
*/

/* Dump canonical translations */
process FETCHPROTEINS {
    input:
        val core_db
    output:
        path "${core_db}_translations.fa", emit: fasta

    beforeScript "ENSCODE=${params.enscode} source ${projectDir}/perl5lib.sh"
    
    script:
    """
    perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname ${core_db} -user ${params.user} -file ${core_db}_translations.fa  ${params.dump_params}
    """

}


/* Create db for linclust */
process CREATEDB {
    input:
        path combined_fasta
    output:
        path "mm_db", emit: linclust_db
    script:
    """
    mkdir -p mm_db
    cd mm_db && mmseqs createdb --dbtype 1 ../$combined_fasta lepidoptera_db
    """
}

process RUNLINCLUST {
    input:
        path linclust_db
    output:
        path "clustered_lepidoptera", emit: clu_lepi
    script: 
    """
    mkdir -p clustered_lepidoptera
    cd clustered_lepidoptera && mmseqs linclust ../$linclust_db/lepidoptera_db clusters temp --min-seq-id 0.95
    """

}

process REPCLUSTER {
    input:
        path clu_lepi
        path linclust_db
    output:
        path "db_cluster_rep", emit: clu_rep
    script:
    """
    mkdir -p db_cluster_rep
    cd db_cluster_rep && mmseqs createsubdb ../$clu_lepi/clusters ../$linclust_db/lepidoptera_db cluster_rep
    """
}

process TOFASTA {
    input:
        path clu_rep
    output:
        path "cluster_rep.fasta", emit: clu_rep_fasta
    script:
    """
    mmseqs convert2fasta $clu_rep/cluster_rep cluster_rep.fasta

    """
}

workflow{
        csvData = Channel.fromPath(params.csvFile).splitCsv().flatten()
        FETCHPROTEINS(csvData)
        combined_fasta = FETCHPROTEINS.out.fasta.collectFile(name: '97_lepidoptera_combined_files.fa', newLine: true)
        CREATEDB(combined_fasta)
        RUNLINCLUST(CREATEDB.out.linclust_db)
        REPCLUSTER(RUNLINCLUST.out.clu_lepi,CREATEDB.out.linclust_db)
        TOFASTA(REPCLUSTER.out.clu_rep)
}
