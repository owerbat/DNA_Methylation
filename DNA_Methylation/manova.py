from statsmodels.multivariate.manova import MANOVA
from pandas import DataFrame
import cpg_sites_to_bops
import numpy as np


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


def correct_cpg_data():
    try:
        annotations = open('..\\methylation_data\\annotations.dat', 'rb')
    except FileNotFoundError:
        print('There is no pickle file. Run create_table script.')
        return 2
    data = cpg_sites_to_bops.pickle.load(annotations)
    data = cpg_sites_to_bops.correct_data(data)
    data = data.sort_values(by='BOP')
    # bop_to_cpg = cpg_sites_to_bops.get_bop_to_cpg_dictionary(data)
    return data


def get_column_dict(cpg_to_bop):
    try:
        file = open('..\\methylation_data\\average_beta.txt', 'r', encoding='utf-8')
    except FileNotFoundError:
        print('There is no age file')
        return 2

    key_errors = 0
    file.readline()
    column_dict = {}
    for line in file:
        line = line.split()
        cpg_name = line.pop(0)
        try:
            bop_name = cpg_to_bop[cpg_name]
            if bop_name not in column_dict:
                column_dict.update({bop_name: [np.array(line, dtype=np.float)]})
            else:
                column_dict[bop_name].append(np.array(line, dtype=np.float))
        except KeyError:
            key_errors += 1
    print('key_errors = ' + str(key_errors))
    file.close()
    return column_dict


def multivariate_anova():
    cpg_data = correct_cpg_data()
    print(1)
    cpg_to_bop = cpg_sites_to_bops.get_cpg_to_bop_dictionary(cpg_data)
    print(2)
    column_dict = get_column_dict(cpg_to_bop)
    print(3)

    del cpg_data
    del cpg_to_bop

    file = open('bop_manova.txt', 'w', encoding='utf-8')
    file.write('BoP_name    p_value\n')
    ages = get_ages()
    p_val_dic = {}
    j = 0
    for bop_name, column_lst in column_dict.items():
        p_val_list = []
        size = len(column_lst)
        if size > 2:
            for i in range(size - 2):
                df = DataFrame({'cpg1': column_lst[i], 'cpg2': column_lst[i + 1], 'cpg3': column_lst[i + 2], 'age': ages})
                model = MANOVA.from_formula('cpg1 + cpg2 + cpg3 ~ age', df)
                test = model.mv_test()
                p_val_list.append(test.results['age']['stat'].values[3, 4])
            minimum = min(p_val_list)
            file.write(bop_name + '\t' + str(minimum) + '\n')
            # p_val_dic.update({bop_name: minimum})
        print(j)
        j += 1

    return p_val_dic


print(multivariate_anova())
