from scipy import stats
import numpy as np


def get_ages():
    file = open("D:\\Projects\\DNA Methylation\\methylation_data\\age.txt", 'r', encoding="utf-8")

    ages = []
    while True:
        current_age = file.readline()
        if len(current_age) == 0:
            break
        ages.append(int(current_age))

    file.close()
    print(ages)
    print(len(ages))
    return ages


def gene_data_read(file):
    gene = file.readline().split(" ")
    gene_name = gene.pop(0)
    print(gene_name)

    for i in range(len(gene)):
        gene[i] = float(gene[i])
    print(gene)
    print(len(gene))

    return gene_name, gene


def main():
    file = open("D:\\Projects\\DNA Methylation\\methylation_data\\gene_data.txt", 'r', encoding="utf-8")

    x = get_ages()
    name, y = gene_data_read(file)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("r_value = " + str(r_value))

    file.close()


main()
