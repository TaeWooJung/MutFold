import mysql.connector
from mysql.connector import Error
# import pandas as pd
import os
from parse_tsv import uniprot_tsv_to_tables, cosmic_tsv_to_tables

# uniprot_tsv_file is from,
# https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue%2Cmodel_organism%3A9606%2Cproteins_with%3A20&query=breast%20cancer

cwd = os.getcwd()
uniprot_tsv_file = 'data/uniprot_breast_cancer.tsv'
cosmic_tsv_file = 'data/cosmic_filtered.tsv'

uniprot_tables = uniprot_tsv_to_tables(os.path.join(cwd, uniprot_tsv_file))
# print(uniprot_tables.keys())
try:
     user = ''
     password = ''

     # connect to local mysql server as a user
     connection = mysql.connector.connect(host='localhost',
                                        user=user,
                                        password=password,
                                        database='cancer_uniprotdb')
     cursor = connection.cursor()

     # remove table if it exists
     cursor.execute("""DROP TABLE IF EXISTS pymol_align""")
     cursor.execute("""DROP TABLE IF EXISTS elaspic""")
     cursor.execute("""DROP TABLE IF EXISTS structure_3d""")
     cursor.execute("""DROP TABLE IF EXISTS active_sites""")
     cursor.execute("""DROP TABLE IF EXISTS binding_sites""")
     cursor.execute("""DROP TABLE IF EXISTS interactions""")
     cursor.execute("""DROP TABLE IF EXISTS keywords""")
     cursor.execute("""DROP TABLE IF EXISTS pubmed""")
     cursor.execute("""DROP TABLE IF EXISTS mutations""")
     cursor.execute("""DROP TABLE IF EXISTS cosmic_data""")
     
     # removing 'proteins' table last since it contains PRIMARY KEY
     cursor.execute("""DROP TABLE IF EXISTS proteins""")

     # create 'proteins' table
     # entry_id, entry_name, protein_names, gene_names, organism, seq, seq_length, dlm
     proteins_colnames, proteins_values = uniprot_tables['proteins']

     cursor.execute("""CREATE TABLE proteins (\
                         entry_id VARCHAR(20) PRIMARY KEY, entry_name VARCHAR(50) NOT NULL, \
                         protein_names VARCHAR(6000) NOT NULL, gene_names VARCHAR(6000) NOT NULL, \
                         organism VARCHAR(50) NOT NULL, seq TEXT, seq_length INT, dlm DATE \
                    )""")
     mysql_insert_query = """INSERT INTO proteins ({}) \
                                   VALUES (%s, %s, %s, %s, %s, %s, %s, %s) """.format(', '.join(proteins_colnames))
     cursor.executemany(mysql_insert_query, proteins_values)

     # create 'active_sites' table
     # entry_id, position, source, evidence_id
     active_sites_colnames, active_sites_values = uniprot_tables['active_sites']
     cursor.execute("""CREATE TABLE active_sites (\
                         entry_id VARCHAR(20), position VARCHAR(20), \
                         evidences VARCHAR(1000), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id) \
                    )""")
     mysql_insert_query = """INSERT INTO active_sites ({}) \
                                   VALUES (%s, %s, %s) """.format(', '.join(active_sites_colnames))
     cursor.executemany(mysql_insert_query, active_sites_values)

     # create 'binding_sites' table
     # entry_id, position, source, evidence_id
     binding_sites_colnames, binding_sites_values = uniprot_tables['binding_sites']
     cursor.execute("""CREATE TABLE binding_sites (\
                         entry_id VARCHAR(20), position VARCHAR(20), ligand NVARCHAR(100), \
                         ligand_id VARCHAR(20), evidences VARCHAR(1000),\
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id) \
                    )""")
     mysql_insert_query = """INSERT INTO binding_sites ({}) \
                                   VALUES (%s, %s, %s, %s, %s) """.format(', '.join(binding_sites_colnames))
     cursor.executemany(mysql_insert_query, binding_sites_values)


     # create 'interactions' table
     # entry_id, interaction
     interactions_colnames, interactions_values = uniprot_tables['interactions']
     cursor.execute("""CREATE TABLE interactions (\
                         entry_id VARCHAR(20), interaction VARCHAR(50), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id) \
                    )""")
     mysql_insert_query = """INSERT INTO interactions ({}) \
                                   VALUES (%s, %s) """.format(', '.join(interactions_colnames))
     
     cursor.executemany(mysql_insert_query, interactions_values)


     # create 'keywords' table
     # entry_id, keyword
     keywords_colnames, keywords_values = uniprot_tables['keywords']
     cursor.execute("""CREATE TABLE keywords (\
                         entry_id VARCHAR(20), keyword VARCHAR(255), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id) \
                    )""")
     mysql_insert_query = """INSERT INTO keywords ({}) \
                                   VALUES (%s, %s) """.format(', '.join(keywords_colnames))
     
     cursor.executemany(mysql_insert_query, keywords_values)


     # create 'pubmed' table
     # entry_id, pubmed_id
     pubmed_colnames, pubmed_values = uniprot_tables['pubmed']
     cursor.execute("""CREATE TABLE pubmed (\
                         entry_id VARCHAR(20), pubmed_id VARCHAR(255), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id) \
                    )""")
     mysql_insert_query = """INSERT INTO pubmed ({}) \
                                   VALUES (%s, %s) """.format(', '.join(pubmed_colnames))

     cursor.executemany(mysql_insert_query, pubmed_values)


     # create 'cosmic_data' table
     # gene_name, mutation_id, mutation_cds, mutation_aa, description 
     cosmic_tables = cosmic_tsv_to_tables(os.path.join(cwd, cosmic_tsv_file), uniprot_tables)
     cosmic_data_colnames, cosmic_data_values = cosmic_tables['cosmic_data']
     
     cursor.execute("""CREATE TABLE cosmic_data (\
                         gene_name VARCHAR(20), \
                         mutation_id VARCHAR(20) PRIMARY KEY, mutation_cds NVARCHAR(6000), \
                         mutation_aa NVARCHAR(6000), description VARCHAR(30) \
                    )""")
     mysql_insert_query = """INSERT INTO cosmic_data ({}) \
                                   VALUES (%s, %s, %s, %s, %s) """.format(', '.join(cosmic_data_colnames))
     
     cursor.executemany(mysql_insert_query, cosmic_data_values)


     # create 'mutations' table
     # gene_name, mutation_id, mutation_cds, mutation_aa, description 
     mutations_colnames, mutations_values = cosmic_tables['mutations']
     
     cursor.execute("""CREATE TABLE mutations (\
                         entry_id VARCHAR(20), \
                         gene_name VARCHAR(20), \
                         mutation_id VARCHAR(20), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id), \
                         FOREIGN KEY (mutation_id) REFERENCES cosmic_data(mutation_id) \
                    )""")
     mysql_insert_query = """INSERT INTO mutations ({}) \
                                   VALUES (%s, %s, %s) """.format(', '.join(mutations_colnames))
     
     cursor.executemany(mysql_insert_query, mutations_values)

     # create 'structure_3d' table
     # entry_id, seq, entry_name, mutation_id, path
     structure_3d_colnames, structure_3d_values = cosmic_tables['structure_3d']
     
     cursor.execute("""CREATE TABLE structure_3d (\
                         entry_id VARCHAR(20), \
                         seq TEXT, \
                         entry_name VARCHAR(50) NOT NULL, \
                         mutation_id VARCHAR(20), \
                         pLDDT FLOAT(24), \
                         pTM FLOAT(24), \
                         path VARCHAR(50), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id) \
                    )""")
     mysql_insert_query = """INSERT INTO structure_3d ({}) \
                                   VALUES (%s, %s, %s, %s, %s, %s, %s) """.format(', '.join(structure_3d_colnames))
     
     cursor.executemany(mysql_insert_query, structure_3d_values)

     # create 'elaspic' table
     # entry_id, mutation_aa, mutation_id, mw_diff, interaction, binding_site, active_site, ddG
     elaspic_colnames, elaspic_values = cosmic_tables['elaspic']
     
     cursor.execute("""CREATE TABLE elaspic (\
                         entry_id VARCHAR(20), \
                         mutation_aa NVARCHAR(50), \
                         mutation_id VARCHAR(20), \
                         mw_diff FLOAT(24), \
                         interaction VARCHAR(20), \
                         binding_site BOOLEAN, \
                         active_site BOOLEAN, \
                         ddG FLOAT(24), \
                         FOREIGN KEY (entry_id) REFERENCES proteins(entry_id), \
                         FOREIGN KEY (mutation_id) REFERENCES cosmic_data(mutation_id) \
                    )""")
     mysql_insert_query = """INSERT INTO elaspic ({}) \
                                   VALUES (%s, %s, %s, %s, %s, %s, %s, %s) """.format(', '.join(elaspic_colnames))
     
     cursor.executemany(mysql_insert_query, elaspic_values)

     # create 'pymol_align' table
     # mutation_id, path
     pymol_align_colnames, pymol_align_values = cosmic_tables['pymol_align']
     
     cursor.execute("""CREATE TABLE pymol_align (\
                         mutation_id VARCHAR(20), \
                         path VARCHAR(50), \
                         FOREIGN KEY (mutation_id) REFERENCES cosmic_data(mutation_id) \
                    )""")
     
     mysql_insert_query = """INSERT INTO pymol_align ({}) \
                                   VALUES (%s, %s) """.format(', '.join(pymol_align_colnames))
     
     cursor.executemany(mysql_insert_query, pymol_align_values)
     connection.commit()

     cursor.execute("""SHOW TABLES""")
     record = cursor.fetchall()
     print(record)

except Error as e:
     print("Error while connecting to MySQL", e)

finally:
     if connection.is_connected():
        cursor.close()
        connection.close()

# get AA residue size from below url
# https://be.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/