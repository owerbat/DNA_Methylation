from scipy import stats
import numpy as np


def get_ages():
    file = open("..\\methylation_data\\age.txt", 'r', encoding="utf-8")

    ages = []
    while True:
        current_age = file.readline()
        if len(current_age) == 0:
            break
        ages.append(int(current_age))

    file.close()

    return ages


def gene_data_read(file):
    gene = file.readline().split(" ")
    gene_name = gene.pop(0)

    for i in range(len(gene)):
        gene[i] = float(gene[i])

    return gene_name, gene


def main():
    file = open("..\\methylation_data\\gene_data.txt", 'r', encoding="utf-8")

    r_values = []
    names = []
    x = get_ages()
    while True:
        name, y = gene_data_read(file)
        if len(y) == 0:
            break
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        r_values.append(r_value)
        names.append(name)

    file.close()

    if len(names) != len(r_values):
        print("Error!")
    else:
        table = list([r_values[i], names[i]] for i in range(len(names)))
        table.sort(reverse=True)

        result_file = open("Results.txt", 'w', encoding="utf-8")
        for i in range(10):
            result_file.write(str(i+1) + ') ' + str(table[i][1]) + ' - ' + str(table[i][0]) + '\n')
        result_file.close()


main()
