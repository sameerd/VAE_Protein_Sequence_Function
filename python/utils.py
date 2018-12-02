def truncate_seq_from_asterisk_onwards(seq):
    """Truncate a protein sequence after we see an asterisk. Also remove the asterisk.
    # Arguments:
        seq (string): string representation of a sequence
    # Returns:
        ret (string): string representation of truncated sequence
    """
    idx = seq.find("*")
    ret = seq
    if idx >= 0:
        ret = seq[:idx]
    return(ret)


def mutate_wt(wt_index, alphabet_sorted, *mutations):
    """Mutate the wild type sequence with mutations. 
    # Arguments:
        wt_index (np.int array): index representation of WT sequence in an ALPHABET
        alphabet_sorted (np.chararray): sorted character array (<U1) of the ALPHABET
        *mutations (string) : Individual mutations. ("A4F" is an example)
    # Returns:
        mutation_index (np.int array): index representation of mutated sequence
    """ 
    mutation_index = wt_index.copy()
    for mutation in mutations:
        replace_char = mutation[-1]
        # convert the one based index of the mutation to the zero based index
        idx = int(mutation[1:-1]) - 1
        mutation_index[idx] = alphabet_sorted.searchsorted(replace_char)
    # FIXME: We do not know how to specifically handle sequences with a * in them. 
    # Lets ignore it for now
    return(mutation_index)
