from scipy import stats


def get_ages():
    try:
        file = open("..\\methylation_data\\age.txt", 'r', encoding="utf-8")
    except FileNotFoundError:
        print("methylation_data.txt is not founded")
        exit(2)

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


def calculate_best_genes(way):
    try:
        file = open("..\\methylation_data\\gene_data.txt", 'r', encoding="utf-8")
    except FileNotFoundError:
        print("gene_data.txt is not founded")
        exit(2)

    r_values = []
    names = []
    x = get_ages()
    while True:
        name, y = gene_data_read(file)
        if len(y) == 0:
            break
        if way == 1:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        elif way == 0:
            r_value, pval = stats.spearmanr(x, y)
        r_values.append(r_value)
        names.append(name)

    file.close()

    if len(names) != len(r_values):
        print("Error!")
    else:
        table = list([r_values[i], names[i]] for i in range(len(names)))
        table.sort(reverse=True)

        if way == 1:
            result_file = open("Results_Linear_Regression.txt", 'w', encoding="utf-8")
        elif way == 0:
            result_file = open("Results_Spearman_Correlation.txt", 'w', encoding="utf-8")
        for i in range(10):
            result_file.write(str(i+1) + ') ' + str(table[i][1]) + ' - ' + str(table[i][0]) + '\n')
        result_file.close()


def main():
    the_way = int(input('''What kind of method do you want to use?
0 - Spearman correlation
1 - Linear regression
... '''))
    if the_way == 0 or the_way == 1:
        calculate_best_genes(the_way)
        print("Done!")
    else:
        print("Wrong value!")


main()
