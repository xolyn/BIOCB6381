

CODON_DICT  = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}


# -------------------------------------------
#       1.1 Reverse strand conversion
# -------------------------------------------
def reverse_complement(seq):
    """
    Args:
        seq (str): a DNA sequence consists of {A,T,C,G}
    Return:
        rev_comp (str): the reverse strand of the given sequence
    """
    comp_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    
    # ---- YOUR CODE STARTS HERE ----
    rev_comp = seq[::-1]
    newseq=""
    for i in range(len(rev_comp)):
        newseq+=comp_dict[rev_comp[i]]
    rev_comp=newseq
    # ----- YOUR CODE ENDS HERE -----

    return rev_comp


# -------------------------------------------
#            1.2 Translation
# -------------------------------------------
def translate(seq: str, codon_dict: dict):
    """
    translate DNA sequence into protein sequence
    Args:
        seq (str): input DNA sequence
        codon_dict (dict): dictionary that maps codons to amino acids
    Returns:
        aa_seq (str): translated protein sequence
    """
    aa_seq = ''
    trseq=""
    # ---- YOUR CODE STARTS HERE ----
    # if seq.startswith("ATG")==True:
    for i in range(0,len(seq),3):
        if seq[i:i+3] in list(codon_dict.keys()):
            trseq+= codon_dict[seq[i:i+3]]

    if "*" in trseq:
        trseq=trseq[:trseq.index("*")]
    
    aa_seq=trseq
    
    # ----- YOUR CODE ENDS HERE -----
    
    return aa_seq


# -------------------------------------------
#           1.3 Global Alignment
# -------------------------------------------
def aln_score_mat(seq1, seq2, score_dict):
    """
    Compute the global alignment score matrix
    Args:
        seq1 (str): sequence-1 please put it vertically in the alignment score matrix
        seq2 (str): sequence-2 please put it horizontally in the alignment score matrix
        score_dict (dict): scoring dictionary with keys {'match', 'mismatch', 'gap'}
    Returns:
        aln_score (list(list)): alignment score matrix (2d-list)
    """
    # ---- Initialize ----
    n_row = len(seq1) + 1 
    n_col = len(seq2) + 1
    aln_score = [[0 for col in range(n_col)] for row in range(n_row)]  # initialize score matrix with 0
    
    # ---- YOUR CODE STARTS HERE ----
    # mat=np.matrix(aln_score)
    for i in range(n_row):
        aln_score[i][0]=i*(score_dict["gap"])
    for i in range(n_col):
        aln_score[0][i]=i*(score_dict["gap"])

    for i in range(0,n_row-1):
        for j in range(0,n_col-1):
            if seq1[i]==seq2[j]:
                aln_score[i+1][j+1]=max(aln_score[i][j]+score_dict["match"],aln_score[i+1][j]+score_dict["gap"],aln_score[i][j+1]+score_dict["gap"])
            else:
                aln_score[i+1][j+1]=max(aln_score[i][j]+score_dict["mismatch"],aln_score[i+1][j]+score_dict["gap"],aln_score[i][j+1]+score_dict["gap"])
    # ----- YOUR CODE ENDS HERE -----
    
    return aln_score


