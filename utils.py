def create_mutated_seq(seq, mutations):

    recognized_mutations = []
    unrecognized_mutations = []
    mutations = mutations.split(';')

    for mutation in mutations:
        
        # trim str
        mutation = mutation[2:]
        aa_before = mutation[0]
        aa_after = mutation[-1]
        pos = int(mutation[1:-1])-1
        
        # whether mutation notation is valid
        if pos < len(seq):

            # if AA matches, change the AA sequence
            if aa_before == seq[pos]:
                recognized_mutations.append('p.{}'.format(mutation))
                seq = seq[:pos] + aa_after + seq[pos+1:]

            else:
                unrecognized_mutations.append('p.{}'.format(mutation))
        
        else:
                unrecognized_mutations.append('p.{}'.format(mutation))

    return seq, recognized_mutations, unrecognized_mutations


def uniprot_gene_to_entry_id_matching(proteins_values):

    gene_to_entry_id = dict()

    for proteins_value in proteins_values:
        
        entry_id = proteins_value[0]
        gene_names = proteins_value[3].split(' ')
        
        for gene_name in gene_names:
            if gene_name in gene_to_entry_id:
                gene_to_entry_id[gene_name].add(entry_id)
            else:
                gene_to_entry_id[gene_name] = {entry_id}

    return gene_to_entry_id


def parse_site_info_to_dict(site_info):

    entry_id_to_site_info = dict()

    for site in site_info:

        entry_id = site[0]
        pos = site[1]

        if not pos.isdigit():
                pos = [int(i) for i in pos.split('..')]
                pos = [i for i in range(pos[0], pos[1]+1)]

        else:
            pos = [int(pos)]

        if entry_id not in entry_id_to_site_info:
                entry_id_to_site_info[entry_id] = pos

        else:
            entry_id_to_site_info[entry_id] += pos

    return entry_id_to_site_info

# seq = 'AGCGG'
# mutations = ['C5C']


def get_difference_in_mw(mutation):
     
    # https://be.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/
    # Amino acid residue molecular weights
    aa_molecular_weights = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147,
                            'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165,
                            'P': 115, 'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117}
    
    aa_before = mutation[0]
    aa_after = mutation[-1]

    return aa_molecular_weights[aa_after] - aa_molecular_weights[aa_before]
     