import pickle
import pandas
import numpy as np
from scipy import stats


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


def prepare__data(to_delete):
    try:
        annotations = open('..\\..\\methylation_data\\annotations.dat', 'rb')
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


def get_ages():
    try:
        file = open('..\\..\\methylation_data\\attributes GSE40279.txt', 'r', encoding='utf-8')
    except FileNotFoundError:
        print('There is no age file')
        return 2

    file.readline()
    ages = []
    male_ages = []
    male_numbers = []
    female_ages = []
    female_numbers = []
    i = 0
    for line in file:
        age = int(line.split()[2])
        ages.append(age)
        sex = line.split()[3]
        if sex == 'M':
            male_ages.append(age)
            male_numbers.append(i)
        else:
            female_ages.append(age)
            female_numbers.append(i)
        i += 1

    file.close()
    return ages, male_ages, female_ages, male_numbers, female_numbers


def get_cpg_data(df, sex=''):
    try:
        data_file = open('..\\..\\methylation_data\\average_beta.txt', 'r', encoding='utf-8')
    except FileNotFoundError:
        print('There is no data file')
        return 2

    cpg_r = {}
    ages, male_ages, female_ages, male_numbers, female_numbers = get_ages()
    data_file.readline()
    j = 0
    for line in data_file:
        y = line.split()
        cpg_name = y.pop(0)
        for i, item in enumerate(y):
            y[i] = float(y[i])
        if sex == '':
            slope, intercept, r_value, p_value, std_err = stats.linregress(ages, y)
        elif sex == 'M':
            male_y = []
            for i in male_numbers:
                male_y.append(y[i])
            slope, intercept, r_value, p_value, std_err = stats.linregress(male_ages, male_y)
        elif sex == 'F':
            female_y = []
            for i in female_numbers:
                female_y.append(y[i])
            slope, intercept, r_value, p_value, std_err = stats.linregress(female_ages, female_y)
        cpg_r.update({cpg_name: r_value})
        print(j)
        j += 1
    print('dictionary is ready')

    j = 0
    r_values = []
    cpg_names = list(df['ID_REF'])
    for cpg in cpg_names:
        try:
            r_values.append(cpg_r[cpg])
        except KeyError:
            j += 1
    print('error count =', j)
    r_values, cpg_names = (list(t) for t in zip(*sorted(zip(r_values, cpg_names), reverse=True)))

    if sex == '':
        cpg_data_file = open('best_cpg.txt', 'w', encoding='utf-8')
    elif sex == 'M':
        cpg_data_file = open('best_male_cpg.txt', 'w', encoding='utf-8')
    elif sex == 'F':
        cpg_data_file = open('best_female_cpg.txt', 'w', encoding='utf-8')
    length = int(len(r_values) * 0.05)
    for i in range(length):
        cpg_data_file.write(cpg_names[i] + '\t' + str(r_values[i]) + '\n')
    cpg_data_file.close()


def cpg_best():
    try:
        best_file = open('best_cpg.txt', 'r', encoding='utf-8')
        male_file = open('best_male_cpg.txt', 'r', encoding='utf-8')
        female_file = open('best_female_cpg.txt', 'r', encoding='utf-8')
    except FileNotFoundError:
        print('There is no data file')
        return 2

    best = set()
    best_dict = {}
    for line in best_file:
        cpg_name = line.split()[0]
        r_value = line.split()[1]
        best.add(cpg_name)
        best_dict.update({cpg_name: r_value})
    best_file.close()

    male = set()
    male_dict = {}
    for line in male_file:
        cpg_name = line.split()[0]
        r_value = line.split()[1]
        male.add(cpg_name)
        male_dict.update({cpg_name: r_value})
    male_file.close()

    female = set()
    female_dict = {}
    for line in female_file:
        cpg_name = line.split()[0]
        r_value = line.split()[1]
        female.add(cpg_name)
        female_dict.update({cpg_name: r_value})
    female_file.close()

    best_male = best & male
    best_female = best & female
    best_best = best_male & best_female
    best_best_male = best_male - best_best
    best_best_female = best_female - best_best

    best_list = []
    for cpg in best_best:
        best_list.append([float(best_dict[cpg]), cpg])
    best_list.sort(reverse=True)
    best_best_file = open('best-best_cpg.txt', 'w', encoding='utf-8')
    for item in best_list:
        best_best_file.write(item[1] + '\t' + str(item[0]) + '\n')
    best_best_file.close()

    best_list = []
    for cpg in best_best_male:
        best_list.append([float(male_dict[cpg]), cpg])
    best_list.sort(reverse=True)
    best_male_file = open('best-best_male_cpg.txt', 'w', encoding='utf-8')
    for item in best_list:
        best_male_file.write(item[1] + '\t' + str(item[0]) + '\n')
    best_male_file.close()

    best_list = []
    for cpg in best_best_female:
        best_list.append([float(female_dict[cpg]), cpg])
    best_list.sort(reverse=True)
    best_female_file = open('best-best_female_cpg.txt', 'w', encoding='utf-8')
    for item in best_list:
        best_female_file.write(item[1] + '\t' + str(item[0]) + '\n')
    best_female_file.close()


def main():
    # data = prepare__data({'UCSC_REFGENE_NAME': [''], 'CHR': ['X', 'Y']})
    # print(data)
    cpg_best()


if __name__ == '__main__':
    main()
