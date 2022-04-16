from itertools import chain

def _get_seq(ranges, x, sep = ' '):
    '''
        A helper function to return a protein subsequence
    '''
    lista = [list( range(valor[0], valor[1] + 1) ) for valor in ranges]
    seq_residues = list(chain.from_iterable(lista))
    seq_residues_str = sep.join(str(e) for e in seq_residues)
    if x == 'str':
        final_seq = seq_residues_str
    elif x == 'list':
        final_seq = seq_residues
    else: 
        final_seq = "Which one? 'str' or 'list'"
    return(final_seq)


def get_pocket_residues(x='str', sep = ' '):
    pocket_rangeResidues = [[8,19], [30,33], [64,65], [79,90], [129,134], [143,146]]
    final_seq = _get_seq(pocket_rangeResidues, x, sep)
    return(final_seq)

def get_secondary_structure_residues(x='str', sep = ' '):
    '''From https://www.uniprot.org/uniprot/P24941'''
    sec_rangeResidues = [[3,12], [14,23], [29,34], [39,41], 
                         [42,44], [46,57], [66,72], [75,81],
                         [84,86], [87,93], [101, 120], 
                         [130, 132], [133,135], [137,139],
                         [141,143], [148,152], [159,161],
                         [166,168], [171, 174], [178,180],
                         [183,198], [208, 219], [230, 232],
                         [248, 251], [252,254], [257,266],
                         [271,273], [277,281], [284,286]]
    final_seq = _get_seq(sec_rangeResidues, x, sep)
    return(final_seq)

def get_pisani_residues(x='str', sep = ' '):
    pisiani_rangeResidues = [ [4,12], [17, 24], [29,34], [46,55], [66,71], [76,81],  
                            [87,93], [101, 120], [121, 135], [140, 150], [182, 194], [277, 282]]
    final_seq = _get_seq(pisiani_rangeResidues, x, sep)
    return(final_seq)