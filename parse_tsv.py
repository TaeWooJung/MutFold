from utils import uniprot_gene_to_entry_id_matching, parse_site_info_to_dict, create_mutated_seq, get_difference_in_mw
import pandas as pd
import os
import re


def uniprot_tsv_to_tables(path):

    # load tsv file into dataframe
    pd.options.display.max_colwidth = None
    df = pd.read_csv(path, sep='\t', header=0)
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
    
    row, _ = df.shape
    tables = dict()
    print('parsing uniprot data from given tsv files...')

    # parsing protein info
    proteins_colnames = ['entry_id', 'entry_name', 'protein_names', 'gene_names', 'organism', 'seq_length', 'seq', 'dlm']
    proteins_values = []
    tmp_df = df.loc[:,['entry_id', 'entry_name', 'protein_names', 'gene_names', 'organism', 'seq_length', 'seq', 'dlm']]
    
    for r in range(row):
        info = list(filter(None, re.split('\n', tmp_df.iloc[r,:].to_string(index=False, header=False))))
        info = [i.strip() for i in info]
        # info[2] = re.split('[\(\)]', info[2])[0]
        info[-3] = int(info[-3])
        proteins_values.append(tuple(info))
        # print(info[0], info[2])
    
    tables['proteins'] = [proteins_colnames, proteins_values]


    # parsing active_sites info
    active_sites_colnames = ['entry_id', 'position', 'evidences']
    active_sites_values = []
    tmp_df = df.loc[:,['active_sites']]

    for r in range(row):
        id = df.iloc[r, 0]
        active_sites = list(filter(None, re.split('ACT_SITE ', tmp_df.iloc[r,:].to_string(index=False, header=False))))
        
        if 'NaN' not in active_sites:
            for active_site in active_sites:
                info = list(filter(None,re.split('; ',active_site)))
                position = info[0]
                evidences = [i for i in info if '/evidence' in i]
                
                for i in range(len(evidences)):
                    evidence = re.search('"(.*?)"', evidences[i]).group(1)
                    output = tuple([id, position, evidence])
                
                active_sites_values.append(output)
        
    tables['active_sites'] = [active_sites_colnames, active_sites_values]


    # parsing binding_sites info
    binding_sites_colnames = ['entry_id', 'position', 'ligand', 'ligand_id', 'evidences']
    binding_sites_values = []
    tmp_df = df.loc[:,['binding_sites']]

    for r in range(row):
        id = df.iloc[r, 0]
        binding_sites = list(filter(None, re.split('BINDING ', tmp_df.iloc[r,:].to_string(index=False, header=False))))
        
        for binding_site in binding_sites:
            position = None
            ligand = None
            ligand_id = None
            evidence = None
    
            if binding_site != 'NaN':
                info = list(filter(None,re.split('; |\|', binding_site)))
                position = info[0]
                ligand = re.search('"(.*?)"', info[1]).group(1)
                try:
                    ligand_id = re.search(':(CHEBI:.*?)"', info[2]).group(1)
                except:
                    ligand_id = None
                # print(info)
                evidences = [i for i in info if '/evidence' in i]
                
                if evidences:
                    evidence_list = []
                    
                    for i in range(len(evidences)):
                        evidence = re.search('/evidence="(.*)', evidences[i]).group(1)
                        evidence_list.append(evidence.strip('"'))

                    evidence = ', '.join(evidence_list)
                
            if position:
                info = tuple([id, position, ligand, ligand_id, evidence])
                binding_sites_values.append(info)
        
    tables['binding_sites'] = [binding_sites_colnames, binding_sites_values]


    # parsing interact_with info
    interactions_colnames = ['entry_id', 'interaction']
    interactions_values = []
    tmp_df = df.loc[:,['interact_with']]

    for r in range(row):
        id = df.iloc[r, 0]
        interactions = list(filter(None, re.split('; ', tmp_df.iloc[r,:].to_string(index=False, header=False))))
        
        for interaction in interactions:
            interaction_id = ''
            if interaction != 'NaN':
                if 'PRO_' in interaction:
                    interaction_id = re.search('\[(.*?)\]', interaction).group(1)
                else:
                    interaction_id = interaction

            if interaction_id:
                info = tuple([id, interaction_id])
                interactions_values.append(info)
    
    tables['interactions'] = [interactions_colnames, interactions_values]

    # parsing keywords info
    keywords_colnames = ['entry_id', 'keyword']
    keywords_values = []
    tmp_df = df.loc[:,['keywords']]

    for r in range(row):
        id = df.iloc[r, 0]
        keywords = list(filter(None, re.split(';', tmp_df.iloc[r,:].to_string(index=False, header=False))))
        
        for keyword in keywords:
            if keyword:
                info = tuple([id, keyword])
                keywords_values.append(info)

    tables['keywords'] = [keywords_colnames, keywords_values]


    # parsing pubmed info
    pubmed_colnames = ['entry_id', 'pubmed_id']
    pubmed_values = []
    tmp_df = df.loc[:,['pubmed_ids']]

    for r in range(row):
        id = df.iloc[r, 0]
        pubmed_ids = list(filter(None, re.split('; ', tmp_df.iloc[r,:].to_string(index=False, header=False))))
        
        for pubmed_id in pubmed_ids:
            if pubmed_id:
                info = tuple([id, pubmed_id])
                pubmed_values.append(info)
    
    tables['pubmed'] = [pubmed_colnames, pubmed_values]
    # print(pubmed_values)

    return tables


def cosmic_tsv_to_tables(path, uniprot_tables, ELASPIC_input=False, ESM_fold_fasta=False):

    # load tsv file into dataframe
    pd.options.display.max_colwidth = None
    df = pd.read_csv(path, sep='\t', header=0)
    # rename headers
    df.rename(columns={"Gene name": "gene_name", 
                    "LEGACY_MUTATION_ID": "mutation_id", 
                    "Mutation CDS": "mutation_cds", "Mutation AA": "mutation_aa", 
                    "Mutation Description": "description"}, inplace=True)
    # print(df.columns)
    print('parsing cosmic data from given tsv files...')
    
    # retrieve gene names availabe from uniprot tables
    uniprot_genes = set()
    _, proteins_values = uniprot_tables['proteins']
    _, binding_sites_values = uniprot_tables['binding_sites']
    _, active_sites_values = uniprot_tables['active_sites']

    for protein_value in proteins_values:
        gene_names = set(protein_value[3].split(' '))
        uniprot_genes = uniprot_genes.union(gene_names)


    # parsing cosmic_mutation info 
    row, _ = df.shape
    cosmic_data_colnames = ['gene_name', 'mutation_id', 'mutation_cds', 'mutation_aa', 'description']
    cosmic_data_values = []
    tables = dict()

    mutation_id = ''
    gene_name = ''
    mutation_cds = []
    mutation_aa = []
    log = []
    output = []

    for r in range(row):
        mutation = list(filter(None, re.split('\n', df.iloc[r,:].to_string(index=False, header=False))))
        mutation = [i.strip() for i in mutation]
        
        # filter cosmic genes if genes do not exist in uniprot
        if mutation[0] in uniprot_genes:
            
            if gene_name != mutation[0]:

                if len(mutation_cds) == 0:
                    mutation_cds.append(mutation[2])
                    mutation_aa.append(mutation[3])
                
                else:
                    info = tuple([gene_name, mutation_id, ';'.join(mutation_cds), ';'.join(mutation_aa), 'Substitution - Missense'])
                    cosmic_data_values.append(info)
                    mutation_cds = [mutation[2]]
                    mutation_aa = [mutation[3]]

                gene_name = mutation[0]
                mutation_id = mutation[1]
            
            else:
                if mutation_id != mutation[1]:
                    info = tuple([gene_name, mutation_id, ';'.join(mutation_cds), ';'.join(mutation_aa), 'Substitution - Missense'])
                    cosmic_data_values.append(info)
                    mutation_cds = [mutation[2]]
                    mutation_aa = [mutation[3]]
                    mutation_id = mutation[1]
                
                else:
                    mutation_cds.append(mutation[2])
                    mutation_aa.append(mutation[3])

    tables['cosmic_data'] = [cosmic_data_colnames, cosmic_data_values]

    # creating mutations table
    gene_to_entry_id = uniprot_gene_to_entry_id_matching(proteins_values)
    mutations_colnames = ['entry_id', 'gene_name', 'mutation_id']
    mutations_values = []
    mutation_id_to_mutation_aa = dict()
    entry_id_to_mutation_id = dict()

    for mutation in cosmic_data_values:
        gene_name = mutation[0]
        mutation_id = mutation[1]
        mutation_aa = mutation[3]
        entry_ids = gene_to_entry_id[gene_name]

        # store mutation info
        mutation_id_to_mutation_aa[mutation_id] = mutation_aa

        for entry_id in entry_ids:
            info = tuple([entry_id, gene_name, mutation_id])
            mutations_values.append(info)
            
            if entry_id not in entry_id_to_mutation_id:
                entry_id_to_mutation_id[entry_id] = {mutation_id}

            else:
                entry_id_to_mutation_id[entry_id].add(mutation_id)

    tables['mutations'] = [mutations_colnames, mutations_values]

    
    entry_id_to_entry_name_and_seq = dict()
    
    for proteins_value in proteins_values:
        entry_id = proteins_value[0]
        entry_name = proteins_value[1]
        seq = proteins_value[6]
        entry_id_to_entry_name_and_seq[entry_id] = [entry_name, seq]
    
    entry_id_to_binding_sites = parse_site_info_to_dict(binding_sites_values)
    entry_id_to_active_sites = parse_site_info_to_dict(active_sites_values)

    entry_id_to_mutation_aa = dict()

    for entry_id in entry_id_to_entry_name_and_seq:
        
        seq = entry_id_to_entry_name_and_seq[entry_id][1]
        entry_name = entry_id_to_entry_name_and_seq[entry_id][0]

        if entry_id in entry_id_to_mutation_id:
            mutation_ids = entry_id_to_mutation_id[entry_id]
            mutation_aa = []

            for mutation_id in mutation_ids:
                mutation_aa = mutation_id_to_mutation_aa[mutation_id]
                
                if entry_id in entry_id_to_mutation_aa:
                    entry_id_to_mutation_aa[entry_id].append([mutation_id, mutation_aa])
                    # entry_name_to_mutation_aa[entry_name] += ';{}'.format(mutation_aa)

                else:
                    entry_id_to_mutation_aa[entry_id] = [[mutation_id, mutation_aa]]
                    # entry_name_to_mutation_aa[entry_name] = mutation_aa
        
    # generated mutated sequence and filter unrecognized mutations
    ELASPIC_entry = []
    ESM_fold_fasta_entry = []

    # parsing 3D protein structure info
    structure_3d_colnames = ['entry_id', 'seq', 'entry_name', 'mutation_id', 'pLDDT', 'pTM', 'path']
    structure_3d_values = []
    
    for entry_id in entry_id_to_mutation_aa:
        seq = entry_id_to_entry_name_and_seq[entry_id][1]
        entry_name = entry_id_to_entry_name_and_seq[entry_id][0]
        mutations = entry_id_to_mutation_aa[entry_id]

        path = 'data/esmfold.log'

        prediction_scores = []

        with open(path, 'r') as f:
            line = f.readline().strip()

            while line:
                line = line.split(',')
                id, pLDDT, pTM  = '', '', ''
                
                for word in line[0].split(' '):
                    if 'HUMAN' in word:
                        id = word + '.pdb'
                if 'pLDDT' in line[1]:
                    pLDDT = line[1].split(' ')[-1]
                    pLDDT = float(pLDDT)
                if 'pTM' in line[2]:
                    pTM = line[2].split(' ')[2]
                    pTM = float(pTM)

                prediction_scores.append((id, pLDDT, pTM))
                line = f.readline().strip()
            f.close()
        
        prediction_scores = pd.DataFrame(prediction_scores, columns=['path', 'pLDDT', 'pTM'])

        # add 3D protein structure info of original protein 
        if len(seq) <= 988:
            prediction_score = prediction_scores.loc[prediction_scores['path'] == '{}.pdb'.format(entry_name),]
            info = tuple([entry_id, seq, entry_name, None, prediction_score['pLDDT'].values[0], prediction_score['pTM'].values[0], '{}.pdb'.format(entry_name)])
            structure_3d_values.append(info)

        if ESM_fold_fasta:
            line = '>{}\n{}'.format(entry_name, seq)
            ESM_fold_fasta_entry.append(line)
        
        for mutation in mutations:
            mutation_id, mutation_aa = mutation
            mutated_seq, recognized_mutations, unrecognized_mutations = create_mutated_seq(seq, mutation_aa)
            # print(mutation_id, ';'.join(recognized_mutations))

            # store valid mutations
            if len(recognized_mutations) != 0:
                for valid_mutation in recognized_mutations:

                    if ELASPIC_input:
                        ELASPIC_entry.append('{}.{}'.format(entry_name, valid_mutation[2:]))

                # was unable to generate 3D structure due to computational issue
                if len(seq) <= 988:
                    structure_3d = '{}_{}.pdb'.format(entry_name, mutation_id)
                    prediction_score = prediction_scores.loc[prediction_scores['path'] == structure_3d,]
                    info = tuple([entry_id, mutated_seq, entry_name, mutation_id, prediction_score['pLDDT'].values[0], prediction_score['pTM'].values[0], structure_3d])
                    structure_3d_values.append(info)

                if ESM_fold_fasta:
                    recognized_mutations = [valid_mutation[2:] for valid_mutation in recognized_mutations]
                    line = '>{}\n{}'.format('{}_{}'.format(entry_name, mutation_id), mutated_seq)
                    ESM_fold_fasta_entry.append(line)

    tables['structure_3d'] = [structure_3d_colnames, structure_3d_values]

    # generate txt file containing entry for ELASPIC
    if ELASPIC_input:
        f = open("data/ELASPIC_entry.txt", "w")
        f.write('\n'.join(ELASPIC_entry))
        f.close()

    # Load ELASPIC result if ELASPIC_input is False
    elif not ELASPIC_input:
        cwd = os.getcwd()
        ELASPIC_path = os.path.join(cwd, 'data/ELASPIC/results.txt')
        pd.options.display.max_colwidth = None
        ELASPIC_df = pd.read_csv(ELASPIC_path, sep='\t', header=0)

        # rename headers
        ELASPIC_df.rename(columns={"Input_identifier": "entry_name", "UniProt_ID": "entry_id", "Mutation": "mutation_aa", 
                        "COSMIC_mut_ID": 'mutation_id', "Interactor_UniProt_ID": "interact_with", "Final_ddG": "ddG",
                        }, inplace=True)
        ELASPIC_df = ELASPIC_df.loc[:,['entry_id', 'entry_name', 'mutation_aa', 'mutation_id', 'interact_with', 'ddG']]

        # parsing ELASPIC info
        elaspic_colnames = ['entry_id', 'mutation_aa', 'mutation_id', 'mw_diff','interaction', 'binding_site', 'active_site', 'ddG']
        elaspic_values = []

        # find interacting entry_id from uniprot 'interactions' table
        interacting_entry_id = []
        for _, interaction in uniprot_tables['interactions'][1]:
            interacting_entry_id.append(interaction)

        for value in ELASPIC_df.values:
            entry_id = value[0]
            mutation_aa = value[2]
            mutation_id = value[3]
            interaction = value[4]

            # if interaction not in interacting_entry_id:
            #     interaction = None
            ddG = value[5]

            # remove information does not have ddG value or mutation_id not in our filtered mutation_id list
            if ddG != '-' and mutation_id in mutation_id_to_mutation_aa:
                if mutation_id == '-':
                    mutation_ids = entry_id_to_mutation_id[entry_id]
                    
                    for x in mutation_ids:
                        if 'p.{}'.format(mutation_aa) in mutation_id_to_mutation_aa[x]:
                            mutation_id = x
                            break

                value = [entry_id, mutation_aa, mutation_id, '', '', '', '', '']
                value = list(map(lambda x: x.replace('-', ''), value))
                # convert ddG to float
                
                value[-5] = get_difference_in_mw(mutation_aa)
                value[-4] = interaction
                value[-3] = False
                value[-2] = False
                value[-1] = float(ddG)

                if entry_id in entry_id_to_binding_sites:
                    if int(mutation_aa[1:-1]) in entry_id_to_binding_sites[entry_id]:
                        value[-3] = True
                    # print(int(mutation_aa[1:-1]), entry_id_to_binding_sites[entry_id])

                if entry_id in entry_id_to_active_sites:
                    if int(mutation_aa[1:-1]) in entry_id_to_active_sites[entry_id]:
                        value[-2] = True
                    # print(int(mutation_aa[1:-1]), entry_id_to_active_sites[entry_id])

                info = tuple(value)
                # print(info)
                elaspic_values.append(info)
            
        tables['elaspic'] = [elaspic_colnames, elaspic_values]

    # write fasta file for ESMfold entry
    if ESM_fold_fasta:
        f = open("data/ESM_fold_entry.fasta", "w")
        f.write('\n'.join(ESM_fold_fasta_entry))
        f.close()

    # parsing pymol_align info
    pymol_align_colnames = ['mutation_id', 'path']
    pymol_align_values = []

    alignments = os.listdir('data/alignments')
    
    mutation_id = ''

    for alignment in alignments:
        if alignment != '.DS_Store':
            mutation_id = alignment.split('_')[-1][:-4]
            pymol_align_values.append(tuple([mutation_id, alignment]))

    tables['pymol_align'] = [pymol_align_colnames, pymol_align_values]
    
    return tables

# cwd = os.getcwd()
# uniprot_tsv_file = 'data/uniprot_breast_cancer.tsv'
# cosmic_tsv_file = 'data/cosmic_filtered.tsv'

# uniprot_tables = uniprot_tsv_to_tables(os.path.join(cwd, uniprot_tsv_file))
# proteins_colnames, proteins_values = uniprot_tables['proteins']
# cosmic_tables = cosmic_tsv_to_tables(os.path.join(cwd, cosmic_tsv_file), uniprot_tables)

