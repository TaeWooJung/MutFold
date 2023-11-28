import pymol
import os
# pymol -cq this_script_name.py.

mutation_id_to_mutation_aa = dict()

with open('./mutation_info.tsv', 'r') as f:
    
    line = f.readline().strip()
    
    while line:
        mutation_id, mutations_aa = line.split('\t')
        mutations_aa = mutations_aa.split(';')
        pos = []
        for mutation_aa in mutations_aa:
            pos.append(mutation_aa[3:-1])
        line = f.readline().strip()
        mutation_id_to_mutation_aa[mutation_id] = pos

    f.close()

cwd = os.getcwd()
structures_path = os.path.join(cwd, 'data/structures')
list_structure_3d = os.listdir(structures_path)

# align_finished = set()

# for png in os.listdir('data/alignments'):
#     png = png.split('_')
    
#     if '.DS' not in png and 'Store' not in png:
#         align_finished.add(f'{png[1]}_{png[2]}')

dict_structure_3d = dict()

for structure in list_structure_3d:
    
    if structure.count('_') > 1:
        gene_name = structure.split('_')
        gene_name = f'{gene_name[0]}_{gene_name[1]}'
        
    else:
        gene_name = structure[:-4]
        
    if gene_name not in dict_structure_3d:
        dict_structure_3d[gene_name] = [structure]

    else:
        dict_structure_3d[gene_name].append(structure)


def super_prots():
    
    path = 'data/structures/'

    for protein in dict_structure_3d:

        prot_1 = f'{path}{protein}.pdb'

        for protein_mut in dict_structure_3d[protein]:

            if '_COSM' in protein_mut:
                prot_2 = f'{path}{protein_mut}'

                mutation_id = prot_2.split('_')[-1][:-4]
                print(prot_2)
                pos = mutation_id_to_mutation_aa[mutation_id]
                pos = '+'.join(pos)

                # If you are using already downloaded pdb files, use load instead:
                pymol.cmd.load(prot_1, 'prot_1')
                pymol.cmd.load(prot_2, 'prot_2')
                pymol.cmd.select(f'SNP', f'id {pos}')

                pymol.cmd.hide('all')
                pymol.cmd.set('ray_opaque_background', 0)

                # Show cartoon and paint them different
                pymol.cmd.show('cartoon', 'prot_1')
                pymol.cmd.color('orange', 'prot_1')
                pymol.cmd.show('cartoon', 'prot_2')
                pymol.cmd.color('lightblue', 'prot_2')
                pymol.cmd.show('sphere', 'SNP')
                pymol.cmd.color('red', 'SNP')

                # Align both proteins
                pymol.cmd.super('prot_1', 'prot_2')

                # Improve image quality, and get png.
                pymol.cmd.set('depth_cue', 0)
                pymol.cmd.set('spec_reflect', 0)
                pymol.cmd.set('cartoon_sampling', 15)
                pymol.cmd.set('ribbon_sampling', 15)
                pymol.cmd.set('antialias', 2)
                pymol.cmd.space('cmyk')

                # This creates high-quality images, but might take a long time to
                # process. Adjust parameters accordingly.

                alignment_name = 'aligned_' + prot_2.split('/')[-1][:-4]
                pymol.cmd.save(f'./data/alignments/{alignment_name}.pdb', state=0)
                print(alignment_name)
                pymol.cmd.png(f'./data/alignments/{alignment_name}.png',
                            width=2000,
                            dpi=300,
                            ray=1
                            )
                pymol.cmd.reinitialize()

if __name__ == '__main__':
    # These are entries from PDB. It can work with already downloaded
    # files, but you have to change the function accordingly
    protein_matrix = [['1EWL', '4W5B'],
                      ['4E2D', '3ATZ'],
                      ]
    for prot_1, prot_2 in protein_matrix:
        super_prots()

