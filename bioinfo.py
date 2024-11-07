import blosum as blo

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
        while (i >= 0 and j >= 0):

    return matrix_score


"""seq1 = input().split()[2]
seq2 = input().split()[2]
gap_pen = int(input().split()[2])"""

for row in calc_score("PTTEINS", "PRTWPSEIN", -1):
    print(" ".join(f"{value[0]:>5}" for value in row))

