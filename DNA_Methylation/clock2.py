import pickle
import pandas
import numpy as np
from scipy import stats
from statsmodels.regression.linear_model import OLS
from sklearn.metrics import mean_absolute_error


pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 1000)


def convert(value):
    try:
        return int(value)
    except ValueError:
        return np.nan


def remove_repeating_elements(lst):
    result = []
    for el in lst:
        if el not in result:
            result.append(el)
    return result


def prepare_data(to_delete):
    try:
        annotations = open('..\\methylation_data\\annotations.dat', 'rb')
    except FileNotFoundError:
        print('There is no pickle file. Run create_table script.')
        return 2

    data = pickle.load(annotations)
    for key, value in to_delete.items():
        for item in value:
            data = data[data[key] != item]

    data['n.CpG'] = data['n.CpG'].apply(convert)
    data.MAPINFO = data.MAPINFO.apply(convert)
    data.UCSC_REFGENE_NAME = data.UCSC_REFGENE_NAME.apply(lambda string: remove_repeating_elements(string.split(';')))
    data.UCSC_REFGENE_ACCESSION = data.UCSC_REFGENE_ACCESSION.apply(lambda string:
                                                                    remove_repeating_elements(string.split(';')))
    data.UCSC_REFGENE_GROUP = data.UCSC_REFGENE_GROUP.apply(lambda string: remove_repeating_elements(string.split(';')))

    annotations.close()

    return data


def get_list_of_cpg(data):
    return list(set(list(data['ID_REF'])))


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


def get_best_cpg(data):
    file = open('..\\methylation_data\\average_beta.txt', 'r', encoding='utf-8')
    file.readline()

    ages = get_ages()
    cpg_r = {}
    j = 0
    for line in file:
        y = line.split()
        cpg_name = y.pop(0)
        y = np.asarray(y, dtype=float)
        slope, intercept, r_value, p_value, std_err = stats.linregress(ages, y)
        cpg_r.update({cpg_name: r_value})
        j += 1
        if j % 100 == 0:
            print(j)

    cpg_list = get_list_of_cpg(data)
    r_values = []
    j = 0
    for cpg in cpg_list:
        try:
            r_values.append([cpg_r[cpg], cpg])
        except KeyError:
            j += 1
    print('errors_count =', j)
    # r_values, cpg_list = (list(t) for t in zip(*sorted(zip(r_values, cpg_list), reverse=True)))
    r_values.sort(reverse=True)

    cpg_data_file = open('best_cpg.txt', 'w', encoding='utf-8')
    length = int(len(r_values) * 0.05)
    for i in range(length):
        # cpg_data_file.write(cpg_list[i] + '\t' + str(r_values[i]) + '\n')
        cpg_data_file.write(str(r_values[i][1]) + '\t' + str(r_values[i][0]) + '\n')
    cpg_data_file.close()


def ordered_beta():
    file = open('best_cpg.txt', 'r', encoding='utf-8')  # cpg
    gene_list = []
    for line in file:
        gene_list.append(line.split()[0])
    file.close()

    file = open('..\\methylation_data\\average_beta.txt', 'r', encoding='utf-8')  # cpg
    tmp = open('tmp1.txt', 'w', encoding='utf-8')
    i = 0
    for line in file:
        if line.split()[0] in gene_list:
            tmp.write(line)
        i += 1
        if i % 100 == 0:
            print(i)
    tmp.close()
    file.close()

    tmp = open('tmp1.txt', 'r', encoding='utf-8')
    gene_dict = {}
    for line in tmp:
        line = line.split()
        name = line.pop(0)
        gene_dict.update({name: line})
    tmp.close()

    file = open('ordered_average_beta.txt', 'w', encoding='utf-8')  # cpg
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
    result_file = open('epigenetic_clock2_errors.txt', 'w', encoding='utf-8')
    global_errors = []
    for gene_count in range(1, 151):
        try:
            file = open('ordered_average_beta.txt', 'r', encoding='utf-8')
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
    data = prepare_data({'UCSC_REFGENE_NAME': [''], 'CHR': ['X', 'Y']})
    epigenetic_clock()


main()
