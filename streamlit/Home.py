import streamlit as st
import mysql.connector
import pandas as pd

from stmol import showmol
import py3Dmol
import biotite.structure.io as bsio
from clipboard import copy

# Page setup
st.set_page_config(page_title="3D Structure Prediction", page_icon="", layout="wide")

st.sidebar.write("##")
st.sidebar.image("streamlit/images/front_page_sb.png")
st.sidebar.write("##")

@st.cache_resource
def init_connection():
    return mysql.connector.connect(**st.secrets["mysql"])

connection = init_connection()

@st.cache_data(ttl=600)
def run_query(query):
    with connection.cursor() as cursor:
        cursor.execute(query)
        return cursor.fetchall()

def get_difference_in_mw(mutation):
     
    # https://be.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/
    # Amino acid residue molecular weights
    aa_molecular_weights = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147,
                            'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165,
                            'P': 115, 'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117}
    
    aa_before = mutation[0]
    aa_after = mutation[-1]

    return aa_molecular_weights[aa_after] - aa_molecular_weights[aa_before]

# stmol
def render_mol(pdb, spin=True, size=(500, 800), zoom=1, background_color='#0e1117'):
    height, width = size
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb,'pdb')
    pdbview.setStyle({'cartoon':{'color':'spectrum'}})
    pdbview.setBackgroundColor(background_color)#('0xeeeeee')
    pdbview.center()
    pdbview.zoomTo()
    pdbview.zoom(zoom, 500)
    pdbview.spin(spin)
    showmol(pdbview, height=height, width=width)

def show_protein(entry_id, df):
    
    st.empty()

    # retrieve protein info db
    query = '''
            SELECT entry_name, seq, seq_length FROM proteins 
            WHERE entry_id = '{}'
            '''.format(entry_id)
    
    entry_name = run_query(query)[0][0]
    seq = run_query(query)[0][1]
    seq_length = run_query(query)[0][2]

    # retrieve protein info from df
    protein_info = df.loc[df['entry_id']==entry_id]

    st.title('[*{0}*](https://www.uniprot.org/uniprotkb/{0}/entry)  |  {1}'.format(entry_id, entry_name))

    # retrieve protein structure info
    query = '''
            SELECT pLDDT, pTM, path FROM structure_3d \
            WHERE entry_id = '{}' AND mutation_id is NULL
            '''.format(entry_id)
    
    pLDDT, pTM, file_name = run_query(query)[0]

    structure_path = 'data/structures/'
    with open(structure_path+file_name) as ifile:
        pdb_string = "".join([x for x in ifile])
        ifile.close()
    
    gene_names = protein_info['gene_names'].values[0]
    st.markdown(f"**Gene name**: {gene_names}")

    protein_names = protein_info['protein_names'].values[0]
    st.markdown(f"**Protein name**: {protein_names}")

    col1, _, col3 = st.columns(3)
    
    with col1:
        render_mol(pdb_string, spin=False, size=(400, 800))
        st.caption(f"*pLDDT = {pLDDT}, **TM-Score = {pTM}")
    
    with col3:
        st.markdown("**Protein Info**:")
        
        with st.expander("**Interact with**"):
            query = '''
                    SELECT interaction FROM interactions \
                    WHERE entry_id = '{}'
                    '''.format(entry_id)
            
            interactions = []

            for row in run_query(query):
                interactions.append(row[0])
            st.markdown(f"**Found {len(interactions)} interacting proteins:**")
            st.markdown(f"{', '.join(interactions)}")
            if len(interactions) > 0:
                st.markdown("Link to [UniProt](https://www.uniprot.org/)")

        with st.expander("**Active Sites**"):
            query = '''
                    SELECT position, evidences FROM active_sites \
                    WHERE entry_id = '{}'
                    '''.format(entry_id)

            active_site_df = pd.DataFrame(run_query(query), columns=['Position', 'Evidences'])
            
            st.markdown(f"**Found {active_site_df.shape[0]} active site(s):**")

            if active_site_df.shape[0] > 0:
                st.dataframe(active_site_df)

        with st.expander("**Binding Sites**"):
            query = '''
                    SELECT position, ligand, evidences FROM binding_sites \
                    WHERE entry_id = '{}'
                    '''.format(entry_id)

            binding_site_df = pd.DataFrame(run_query(query), columns=['Position', 'Ligand', 'Evidences'])
            
            st.markdown(f"**Found {binding_site_df.shape[0]} binding site(s):**")

            if binding_site_df.shape[0] > 0:
                st.dataframe(binding_site_df)
        
        with st.expander("**PubMed ID**"):
            query = '''
                    SELECT pubmed_id FROM pubmed \
                    WHERE entry_id = '{}'
                    '''.format(entry_id)
            
            pubmed_ids = []

            for row in run_query(query):
                pubmed_ids.append(row[0])
            st.markdown(f"**Found {len(pubmed_ids)} papers related to the protein:**")
            st.markdown(f"{', '.join(pubmed_ids)}")
            st.markdown("Link to [PubMed](https://pubmed.ncbi.nlm.nih.gov/)")
        
        with st.expander("**Keywords**"):
            keywords = ', '.join(protein_info['keywords'].values[0].split(';'))
            st.markdown(f"{keywords}")
    
    st.caption('*pLDDT: is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
    st.caption('**TM-Score: is a metric for assessing the topological similarity of protein structures on a scale from 0-1.')

    st.markdown(f'Protein Sequence ({seq_length}bp):')
    st.code(seq, language='python')

    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )

    # retrieve mutation info
    query = '''
            SELECT mutation_id FROM structure_3d \
            WHERE entry_id = '{}' AND mutation_id is NOT NULL
            '''.format(entry_id)

    mutation_ids = run_query(query)

    st.header("Mutation Info:")

    active_sites = list(map(int, (active_site_df["Position"].values)))
    
    binding_sites_temp = list(binding_site_df["Position"].values)
    binding_sites = []

    for pos in binding_sites_temp:
        if '..' in pos:
            lower, upper = pos.split('..')
            binding_sites += [i for i in range(int(lower), int(upper)+1)]

        else:
            binding_sites.append(int(pos))

    mutation_info_df = show_mutation_summary(entry_id, entry_name, active_sites, binding_sites)

    st.caption(f"**{len(mutation_ids)} mutation(s) found:**")

    for mutation_id in mutation_ids:
        show_mutation(mutation_id[0], mutation_info_df)


def show_mutation_summary(entry_id, entry_name, active_sites, binding_sites):
    
    query = '''
            SELECT mutation_aa, mutation_id FROM cosmic_data \
            WHERE mutation_id IN
            (SELECT mutation_id FROM mutations WHERE entry_id = '{}')
            '''.format(entry_id)
    
    active_site, binding_site = 'X', 'X'
    mutation_info = []
    mutation_types = dict()
    # mutation_summary = dict()
    # i = 0
    for row in run_query(query):
        mutations = row[0].split(';')
        mutation_id = row[1]
        
        for mutation in mutations:
            
            mutation_type = f'{mutation[2]} > {mutation[-1]}'

            if mutation_type in mutation_types:
                mutation_types[mutation_type] += 1
            else:
                mutation_types[mutation_type] = 1

            mw_diff = get_difference_in_mw(mutation[2:])
            pos = int(mutation[3:-1])

            if pos in active_sites:
                active_site = 'O'
            else:
                active_site = 'X'
            
            if pos in binding_sites:
                binding_site = 'O'
            else:
                binding_site = 'X'
            
            mutation_info.append([mutation_id, mutation[2:], mw_diff, active_site, binding_site])

    mutation_info_df = pd.DataFrame(mutation_info, columns=['Mutation ID', 'SNPs', 'MW difference', 'Active site', 'Binding site'])
    
    st.subheader(f"**Mutation summary of {entry_name}:**")
    mutation_types_df = pd.DataFrame.from_dict(mutation_types, orient='index', columns=['Frequency'])

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Frequency of different SNPs in the protein**")
        st.bar_chart(mutation_types_df, y='Frequency', use_container_width=False)
    
    with col2:
        st.markdown("**Characteristics of different SNPs**")
        st.dataframe(mutation_info_df, width=600)
        st.caption("**Note:** Active site & Binding site columns denote whether certain SNP is located in either in active or binding sites")
    
    return mutation_info_df


def show_mutation(mutation_id, mutation_info_df):
    
    alignment_path = 'data/alignments/'

    query = '''
            SELECT mutation_aa FROM cosmic_data \
            WHERE mutation_id = '{}'
            '''.format(mutation_id)
    
    mutation_aa_str = run_query(query)[0][0]
    mutation_aa_str = mutation_aa_str.replace(';', ', ')

    with st.expander(f'**{mutation_id}_{mutation_aa_str}**'):

        # retrieve protein structure info
        query = '''
                SELECT pLDDT, pTM, path FROM structure_3d \
                WHERE mutation_id = '{}'
                '''.format(mutation_id)
        
        pLDDT, pTM, file_name = run_query(query)[0]

        structure_path = 'data/structures/'
        print(structure_path)
        with open(structure_path+file_name) as ifile:
            pdb_string = "".join([x for x in ifile])
            ifile.close()
        
        st.header(f"SNPs: {mutation_aa_str}")

        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader(f"**3D Structure**:")
            render_mol(pdb_string, spin=False, size=(500, 800), zoom=1)
            st.caption(f"pLDDT = {pLDDT}, TM-Score = {pTM}")

        with col2:
            st.subheader(f"**PyMol Alignment:**")
            query= '''
                    SELECT path FROM pymol_align \
                    WHERE mutation_id = '{}'
                    '''.format(mutation_id)
            
            img_name = run_query(query).pop()[0]
            # caption = 'Orange: Wild Type, Lightblue: Mutant, Red: SNPs'
            st.image(alignment_path+img_name, width=600)
            legend_file = 'streamlit/images/alignment_legend.png'
            st.image(legend_file, width=200)
        
        col1, col2 = st.columns(2)

        with col1: 
            st.download_button(
                    label="Download PDB",
                    data=pdb_string,
                    file_name=file_name,
                    mime='text/plain',
                )
            
        with col2:
            st.download_button(
                    label="Download PNG",
                    data=open(alignment_path+img_name, 'rb').read(),
                    file_name=img_name,
                    mime='text/plain',
                )
        
        st.subheader(f"**Characteristics of SNPs:**")
        st.dataframe(mutation_info_df.loc[mutation_info_df['Mutation ID']==mutation_id,], width=600)
        st.caption("**Note:** Active site & Binding site columns denote whether certain SNP is located in either in active or binding sites")

        query= '''
                SELECT mutation_id, mutation_aa, interaction, mw_diff, ddG FROM elaspic \
                WHERE mutation_id = '{}'
                '''.format(mutation_id)
        
        output = run_query(query)

        if len(output) > 0:
            st.subheader(f"**ELASPIC results:**")
            elaspic_df = pd.DataFrame(output, columns=['Mutation ID', 'SNPs', 'Interaction', 'MW difference', '*ddG (kcal/mol)'])
        
            st.dataframe(elaspic_df, width=600)

            st.caption("*ddG: The final Gibbs free energy change in kcal/mol is predicted using more than 70 sequential, molecular and\
                        energentic features with the Stochastic Gradient Boosting of Decision Trees algorihm.")
            st.markdown("[UniProt search](https://www.uniprot.org/)")
            st.caption("Visit [ELASPIC website](http://elaspic.kimlab.org/help/) for more information.")


# Load uniprot entry_ids which can be visualized
query = '''
        SELECT entry_id FROM proteins 
        WHERE entry_id IN 
        (SELECT entry_id FROM structure_3d WHERE mutation_id is NULL)
        '''

updated_entry_id = [i[0] for i in run_query(query)]

# Load tsv file into dataframe
pd.options.display.max_colwidth = None
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
text_search = st.sidebar.text_input(label="Search proteins by UniProt entry id or keywords (e.g., O96017)", value="O96017", )

# Filter the dataframe using masks
entry_id = df_temp["entry_id"].str.contains(text_search)
entry_name = df_temp["entry_name"].str.contains(text_search)
protein_names = df_temp["protein_names"].str.contains(text_search)
gene_names = df_temp["gene_names"].str.contains(text_search)
interact_with = df_temp["interact_with"].str.contains(text_search)
keywords = df_temp["keywords"].str.contains(text_search)
pubmed_ids = df_temp["pubmed_ids"].str.contains(text_search)

df_search = df_temp[entry_id | entry_name | protein_names | gene_names | interact_with | keywords | pubmed_ids]

# count number of mutations for each entry_id
entry_ids = df_search.iloc[:,0]
mutation_n = []

for entry_id in entry_ids:
    # Count number of mutations available for the entry_id
    query = '''
            SELECT COUNT(mutation_id) FROM structure_3d \
            WHERE mutation_id IS NOT NULL AND entry_id = '{}' \
            '''.format(entry_id)
    mutation_n.append(run_query(query).pop()[0])

df_search = df_search.assign(mutation_n=mutation_n)
df_search.sort_values(by=['mutation_n'], ascending=False, inplace=True)
# Show the results, if you have a text_search
# Another way to show the filtered results
# Show the cards
N_cards_per_row = 1
buttons = []
entry_ids = []

# if 'button' not in st.session_state:
#     st.session_state['button'] = 0

if text_search:
    st.sidebar.caption(f"{df_search.shape[0]} result(s) found.")
    st.session_state['button'] = 0

    for n_row, row in df_search.reset_index().iterrows():
        i = n_row%N_cards_per_row
        if i==0:
            cols = st.columns(N_cards_per_row)
            st.sidebar.write("---")
            
        # draw the card
        with cols[n_row%N_cards_per_row]:
            st.sidebar.caption(f"Last modified: {row['dlm'].strip()}")
            st.sidebar.markdown(f"{row['entry_id'].strip()} | {row['entry_name'].strip()} | No. mutation(s): {row['mutation_n']}")
            buttons.append(st.sidebar.button('Show protein', key=row['entry_id'].strip()))
            entry_ids.append(row['entry_id'].strip())
            # st.caption(f"{row['entry_name'].strip()} | No. mutations: {mutation_n}")
            # st.markdown(f"**{row['Video']}**")
    
    for i, b in enumerate(buttons): 
        if buttons[i]:
            st.session_state['button'] = i+1

    if st.session_state['button'] > 0:
        entry_id =  entry_ids[st.session_state['button']-1]
        print(entry_id)
        # st.write(f"Button {st.session_state['button']} was clicked.")
        show_protein(entry_id, df_search)


else:
    st.warning('👈 Search protein(s) to visualize!')


# Show the dataframe (we'll delete this later)
# st.write(df_search)


