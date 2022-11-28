import sqlalchemy
from sqlalchemy import create_engine, MetaData, Table
import pandas as pd
import requests
#import wget
#print(sqlalchemy.__version__)

df = pd.read_excel("Lepidoptera_protein_files.xlsx")

def get_dburl(sql_server_details, core_db_name):
    dburl = sql_server_details + core_db_name
    return dburl

def construct_ftp_url(species, assembly, date, pep_file):
    ftp_url = "https://ftp.ensembl.org/pub/rapid-release/species/" + species + "/" + assembly + "/" + "geneset/" + date + "/" + pep_file
    return ftp_url

def construct_pep_file_name(info_tuple):
    return "-".join(info_tuple) + "-" + "pep.fa.gz"

for index, row in df.iterrows():

    sps = "_".join(row["Species"].split(" "))
    #print(sps, row["GCA"], row["Core db on mysql-ens-sta-5"])
    dburl = get_dburl("mysql://ensro@mysql-ens-sta-5:4684/",row["Core db on mysql-ens-sta-5"])
    print(dburl)
    engine = create_engine(dburl,echo = True)
    conn = engine.connect()
    stmt = 'SELECT meta_value FROM meta WHERE meta_key ="genebuild.last_geneset_update"'
    result_proxy = conn.execute(stmt)
    results = result_proxy.fetchall()
    date = "_".join((results[0][0]).split("-"))
    print(date)
    pep_file = construct_pep_file_name((sps,row["GCA"],date))
    print(pep_file)
    ftp_link = construct_ftp_url(sps,row["GCA"],date,pep_file )
    #ftp_link = construct_ftp_url(sps,row["GCA"],date)
    print(ftp_link)
    response = requests.get(ftp_link)
    with open(pep_file, 'wb') as fd:
        for chunk in response.iter_content(chunk_size=128):
            fd.write(chunk)
      

    




#print(engine.table_names())
#metadata = MetaData()
#meta_table = Table('meta', metadata, autoload=True, autoload_with=engine)
#print(repr(meta_table))

    
# stmt_url = 'SELECT meta_value FROM meta WHERE meta_key ="species.url"'
# result_proxy1 = conn.execute(stmt_url)
# results_url = result_proxy1.fetchall()
# url = results_url

