import blosum as blo
from Bio import pairwise2
#from Bio.SubsMat.MatrixInfo import blosum62
global matrix
matrix = blo.BLOSUM(62)


def calc_score(seq1, seq2, gapp):
    matrix_score = [[]]
    matrix_score[0].append([0, "S"])  # START
    for i in range(1, len(seq1) + 1):
        matrix_score[0].append([matrix_score[0][i - 1][0] + gapp, "I"])  # INSERTION
    for j in range(1, len(seq2) + 1):
        matrix_score.append([[matrix_score[j - 1][0][0] + gapp, "D"]])  # DELETION
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            match = matrix_score[i - 1][j - 1][0] + matrix[seq1[j - 1]][seq2[i - 1]]
            delet = matrix_score[i - 1][j][0] + gapp
            inst = matrix_score[i][j - 1][0] + gapp
            values = [match, delet, inst]
            max_value = max(values)
            index_of_max = values.index(max_value)
            where = "blank"
            if index_of_max == 0:
                where = "M"
            elif index_of_max == 1:
                where = "D"
            elif index_of_max == 2:
                where = "I"
            matrix_score[i].append([max_value, where])
    target = ""
    alignment = ""
    query = ""
    i = len(seq2)
    j = len(seq1)
    while (i >= 1 and j >= 1):
        if matrix_score[i][j][1] == "I":
            target = seq1[j-1] + target
            alignment = "-" + alignment
            query = "-" + query
            j-=1
        elif matrix_score[i][j][1] == "D":
            target = "-" + target
            alignment = "-" + alignment
            query = seq2[i-1] + query
            i-=1
        else:
            target = seq1[j-1] + target
            if seq1[j-1] == seq2[i-1]:
                alignment = "|" + alignment
            else:
                alignment = "." + alignment
            query = seq2[i-1] + query
            j-=1
            i-=1
    return matrix_score[len(seq2)][len(seq1)][0], "".join("target  0", target, len(seq1), "        0", len(seq2), alignment, len(alignment), "query   0",query, len(seq2))


"""seq1 = input().split()[2]
seq2 = input().split()[2]
gap_pen = int(input().split()[2])

for row in calc_score("PTTEINS", "PRTWPSEIN", -1):
    print(" ".join(f"{value[0]:>5}" for value in row))"""
seq1 = "PTTEINS"
seq2 = "PRTWPSEIN"
score, target, alignment, query = calc_score(seq1, seq2, -1)
print(score)
print("target  0",target, len(seq1))
print("        0", alignment, len(alignment))
print("query   0",query, len(seq2))

#alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -1, -1)
#print(alignments)

