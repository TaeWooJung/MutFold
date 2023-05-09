import streamlit as st
import mysql.connector
import pandas as pd
from navigate import nav_page

# Page setup
st.set_page_config(page_title="3D Structure Prediction", page_icon="", layout="wide")
st.image("streamlit/images/front_page.png")
st.write("##")
st.write("##")

@st.cache_resource
def init_connection():
    return mysql.connector.connect(**st.secrets["mysql"])

connection = init_connection()

@st.cache_data(ttl=600)
def run_query(query):
    with connection.cursor() as cursor:
        cursor.execute(query)
        return cursor.fetchall()

# Load uniprot entry_ids which can be visualized
query = '''
        SELECT entry_id FROM proteins \
        WHERE entry_id IN \
        (SELECT entry_id FROM structure_3d WHERE mutation_id is NULL)'''

updated_entry_id = [i[0] for i in run_query(query)]

# Load tsv file into dataframe
uniprot_tsv_file = "data/uniprot_breast_cancer.tsv"
df = pd.read_csv(uniprot_tsv_file, sep='\t', header=0)

# rename headers
df.rename(columns={"Entry": "entry_id", "Entry Name": "entry_name", "Protein names": "protein_names", 
                "Length": "seq_length", "Organism": "organism", "Date of last modification": "dlm",
                "Gene Names": "gene_names",
                "Active site": "active_sites", 
                "Binding site": "binding_sites", 
                "Interacts with": "interact_with",
                "Keywords": "keywords",
                "PubMed ID": "pubmed_ids",
                "Sequence": "seq"
                }, inplace=True)

df_temp = df[df["entry_id"].isin(updated_entry_id)]
df_temp = df_temp.loc[:,['entry_id', 'entry_name', 'protein_names', 'gene_names', 'interact_with', 'keywords', 'pubmed_ids', 'dlm']]

# Use a text_input to get the keywords to filter the dataframe
text_search = st.text_input(label="Search proteins by UniProt entry id or keywords (e.g., Q6P4A7)", value="Q6P4A7", )

# Filter the dataframe using masks
entry_id = df_temp["entry_id"].str.contains(text_search)
entry_name = df_temp["entry_name"].str.contains(text_search)
protein_names = df_temp["protein_names"].str.contains(text_search)
gene_names = df_temp["gene_names"].str.contains(text_search)
interact_with = df_temp["interact_with"].str.contains(text_search)
keywords = df_temp["keywords"].str.contains(text_search)
pubmed_ids = df_temp["pubmed_ids"].str.contains(text_search)

df_search = df_temp[entry_id | entry_name | protein_names | gene_names | interact_with | keywords | pubmed_ids]



entry_ids = df_search.iloc[:,0]
mutation_n = []

# count
for entry_id in entry_ids:
    print(entry_id)
    # Count number of mutations available for the entry_id
    query = '''
            SELECT COUNT(mutation_id) FROM structure_3d \
            WHERE mutation_id IS NOT NULL AND entry_id = '{}' \
            '''.format(entry_id)
    mutation_n.append(run_query(query).pop()[0])

df_search['mutation_n'] = mutation_n
df_search.sort_values(by=['mutation_n'], ascending=False, inplace=True)
# Show the results, if you have a text_search
# Another way to show the filtered results
# Show the cards
N_cards_per_row = 1

if text_search:
    st.caption(f"{df_search.shape[0]} result(s) found.")
    st.write("---")

    for n_row, row in df_search.reset_index().iterrows():
        i = n_row%N_cards_per_row
        if i==0:
            cols = st.columns(N_cards_per_row)
            st.write("---")
            
        # Count number of mutations available for the entry_id
        # query = '''
        #         SELECT COUNT(mutation_id) FROM structure_3d \
        #         WHERE mutation_id IS NOT NULL AND entry_id = '{}' \
        #         '''.format(row['entry_id'])
        # mutation_n = run_query(query).pop()[0]

        # draw the card
        with cols[n_row%N_cards_per_row]:
            st.caption(f"Last modified: {row['dlm'].strip()}")
            st.markdown(f"{row['entry_id'].strip()} | {row['entry_name'].strip()} | No. mutation(s): {row['mutation_n']}")
            if st.button('Show protein', key=row['entry_id'].strip()):
                nav_page("Protein")
            # st.caption(f"{row['entry_name'].strip()} | No. mutations: {mutation_n}")
            # st.markdown(f"**{row['Video']}**")

# Show the dataframe (we'll delete this later)
# st.write(df_search)


