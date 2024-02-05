# import pandas
# import numpy
# import os
# script_dir = os.path.dirname(__file__)
# os.chdir(script_dir)

# # path=os.path.join(os.path.join(".","data"),"result.csv")
# a=pandas.read_csv("result.csv")

# print(a[a["Master Protein Accessions"]=="Q06587"].count()[0])
import numpy as np
seq1 = 'GATTACA'
seq2 = 'GTCGACGCA'
score_dict = {'match': 1, 'mismatch': -1, 'gap': -2}
n_row = len(seq1) + 1 
n_col = len(seq2) + 1
aln_score = [[0 for col in range(n_col)] for row in range(n_row)]  # initialize score matrix with 0
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

print(np.matrix(aln_score))
