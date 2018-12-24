from matplotlib import pyplot as plt
import numpy as np
from statsmodels.regression.linear_model import OLS
from sklearn.metrics import mean_absolute_error


def get_ages():
    try:
        file = open('..\\methylation_data\\attributes GSE40279.txt', 'r', encoding='utf-8')
    except FileNotFoundError:
        print('There is no age file')
        return 2

    file.readline()
    ages = []
    for line in file:
        ages.append(int(line.split()[2]))

    file.close()
    return ages


def get_age_histogram(file_name):
    file = open(file_name, 'r', encoding='utf-8')
    if file_name == '..\\methylation_data\\attributes GSE40279.txt':
        age_column_number = 2
        sex_column_number = 3
        title = 'GSE40279'
        figure_number = 1
    elif file_name == '..\\methylation_data\\attributes GSE87571.txt':
        age_column_number = 3
        sex_column_number = 2
        title = 'GSE87571'
        figure_number = 2
    else:
        print('File error')
        exit(1)

    ages = []
    male_ages = []
    female_ages = []
    file.readline()
    for line in file:
        age = int(line.split()[age_column_number])
        ages.append(age)

        sex = line.split()[sex_column_number]
        if sex == 'M':
            male_ages.append(age)
        else:
            female_ages.append(age)

    plt.figure(figure_number)
    plt.xlabel('Age')
    plt.ylabel('Count')
    plt.title(title)
    plt.hist(ages, bins=np.arange(min(ages)+0.5, max(ages), 1), density=True, histtype='step', color='black')
    plt.hist(male_ages, bins=np.arange(min(ages) + 0.5, max(ages), 1), density=True, histtype='step', color='blue')
    plt.hist(female_ages, bins=np.arange(min(ages) + 0.5, max(ages), 1), density=True, histtype='step', color='pink')
    plt.show()


def ordered_beta():
    # file = open('Epigenetic clock\\gene\\best-best_gene.txt', 'r', encoding='utf-8')  # gene
    file = open('Epigenetic clock\\cpg\\best_cpg.txt', 'r', encoding='utf-8')  # cpg
    gene_list = []
    for line in file:
        gene_list.append(line.split()[0])
    file.close()

    # file = open('Results\\gene_average_beta.txt', 'r', encoding='utf-8')  # gene
    file = open('..\\methylation_data\\average_beta.txt', 'r', encoding='utf-8')  # cpg
    tmp = open('tmp.txt', 'w', encoding='utf-8')
    i = 0
    for line in file:
        if line.split()[0] in gene_list:
            tmp.write(line)
        i += 1
        if i % 100 == 0:
            print(i)
    tmp.close()
    file.close()

    tmp = open('tmp.txt', 'r', encoding='utf-8')
    gene_dict = {}
    for line in tmp:
        line = line.split()
        name = line.pop(0)
        gene_dict.update({name: line})
    tmp.close()

    # file = open('ordered_gene_average_beta.txt', 'w', encoding='utf-8')  # gene
    file = open('ordered_cpg_average_beta_best.txt', 'w', encoding='utf-8')  # cpg
    i = 0
    for gene in gene_list:
        file.write(gene)
        try:
            for item in gene_dict[gene]:
                file.write(' ' + item)
            file.write('\n')
        except KeyError:
            i += 1
    print('key errors:', i)
    file.close()


def epigenetic_clock():
    result_file = open('Results\\epigenetic_clock_errors.txt', 'w', encoding='utf-8')
    global_errors = []
    for gene_count in range(1, 151):
        try:
            # file = open('Results\\ordered_gene_average_beta.txt', 'r', encoding='utf-8')  # gene
            file = open('Results\\ordered_cpg_average_beta_best.txt', 'r', encoding='utf-8')
        except FileNotFoundError:
            print('File not founded')
            return 2

        ages = get_ages()
        signs = np.zeros(shape=(gene_count, 656))
        for i in range(gene_count):
            line = file.readline().split()
            signs[i] = np.asarray(line[1:], dtype=float)
        signs = signs.transpose()
        iter_count = 4
        size = int(656 / iter_count)

        errors = []
        for i in range(iter_count):
            test = np.asarray(signs[size*i: size*(i+1)])
            x = np.asarray(list(signs[:size*i]) + list(signs[size*(i+1):]))
            y = np.asarray(list(ages[:size*i]) + list(ages[size*(i+1):]))

            model = OLS(y, x)
            result = model.fit()
            predicted_ages = result.predict(test)
            errors.append(mean_absolute_error(ages[size*i: size*(i+1)], predicted_ages))
        err = sum(errors)/iter_count
        print('Error' + str(gene_count) + ':', err)
        result_file.write(str(gene_count) + '\t' + str(err) + '\n')
        global_errors.append(err)

        file.close()

    numbers = range(1, 151)
    global_errors, numbers = (list(t) for t in zip(*sorted(zip(global_errors, numbers))))
    print('Best:')
    for i in range(10):
        print(str(numbers[i]) + '\t' + str(global_errors[i]))

    result_file.close()


def main():
    # get_age_histogram('..\\methylation_data\\attributes GSE40279.txt')
    # get_age_histogram('..\\methylation_data\\attributes GSE87571.txt')
    epigenetic_clock()


main()
