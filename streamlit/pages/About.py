import streamlit as st

# Page setup
st.set_page_config(page_title="3D Structure Prediction", page_icon="", layout="wide")

st.title('About MutFold')
st.subheader("**Purpose and Tools**")
st.write('**MutFold** was created to visualize & compare the 3D protein structure of mutant proteins related to human breast cancer. \
         **Evolutionary Scale Modeling (ESMFold)** from Meta ([Github page](https://github.com/facebookresearch/esm)) was used to predict 3D protein \
         structures for both wild types and mutants. Alignments between mutant and wild-type proteins were performed using [**PyMol**](https://pymol.org/2/) software.\
         To further assess the effect of mutation on protein, [**ELASPIC**](http://elaspic.kimlab.org/help/) tool was used to predict the impact of protein affinity towards related proteins.')

st.subheader("**Resoruces**")
st.write('MutFold contains a 3D protein structure of 2314 proteins (251 wild-type proteins (protein length < 989bp) and 2063 mutant proteins. \
          Proteins with sequences longer than 989bp showed a rapid decrease in prediction score, hence, limiting the length to 988bp. However, it is important to consider that \
          the prediction was done using the default parameters of the model provided by the ESMFold repository and parameters might not be optimal for longer sequences. \
          All proteins related to breast cancer and mutations were retrieved from [**UniProt**](https://www.uniprot.org/) and [**COSMIC**](https://cancer.sanger.ac.uk/cosmic) respectively.\
          For the database, MySQL was used to store information gathered from above resources.')

st.subheader("**Schema**")
schema_path = 'streamlit/images/schema.png'
st.image(schema_path, width=700)

st.subheader("**Download Files**")
st.caption("**UniProt:**")
uniprot_file = 'data/uniprot_breast_cancer.tsv'
st.download_button(
            label="Download File",
            data=open(uniprot_file, 'rb').read(),
            file_name='uniprot_breast_cancer.tsv',
            mime='text/plain',
        )
st.caption("**COSMIC:**")
cosmic_file = 'data/cosmic_filtered.tsv'
st.download_button(
            label="Download File",
            data=open(cosmic_file, 'rb').read(),
            file_name='cosmic_filtered.tsv',
            mime='text/plain',
        )

st.caption("**ELASPIC:**")
elaspic_entry_file = 'data/ELASPIC_entry.txt'
elaspic_file = 'data/ELASPIC/results.txt'
st.download_button(
            label="Download Entry File",
            data=open(elaspic_entry_file, 'rb').read(),
            file_name='ELASPIC_entry.txt',
            mime='text/plain',
        )
st.download_button(
            label="Download Result File",
            data=open(elaspic_file, 'rb').read(),
            file_name='ELASPIC_results.txt',
            mime='text/plain',
        )
    
st.caption("**ESMFold:**")
fasta_file = 'data/ESM_fold_entry_filtered_988.fasta'
st.download_button(
            label="Download Fasta File",
            data=open(fasta_file, 'rb').read(),
            file_name='ESM_fold_entry_filtered_988.fasta',
            mime='text/plain',
        )
