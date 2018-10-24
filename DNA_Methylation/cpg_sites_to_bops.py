import pandas
import pickle
import numpy as np


pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 1000)
# pandas.set_option('display.max_rows', 300)


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


def correct_data(data):
    names = list(data.columns)
    names.pop(3)
    names.pop(3)

    data = data[names]
    data = data[data.UCSC_REFGENE_NAME != '']
    data = data[data.RELATION_TO_UCSC_CPG_ISLAND != 'N_Shelf']
    data = data[data.RELATION_TO_UCSC_CPG_ISLAND != 'S_Shelf']
    data = data[data.RELATION_TO_UCSC_CPG_ISLAND != '']
    data = data[data.CHR != 'X']
    data = data[data.CHR != 'Y']
    data = data[data.CHR != '']

    data['n.CpG'] = data['n.CpG'].apply(convert)
    data.MAPINFO = data.MAPINFO.apply(convert)
    data.UCSC_REFGENE_NAME = data.UCSC_REFGENE_NAME.apply(lambda string: remove_repeating_elements(string.split(';')))
    data.UCSC_REFGENE_ACCESSION = data.UCSC_REFGENE_ACCESSION.apply(lambda string:
                                                                    remove_repeating_elements(string.split(';')))
    data.UCSC_REFGENE_GROUP = data.UCSC_REFGENE_GROUP.apply(lambda string: remove_repeating_elements(string.split(';')))

    return data


def get_cpg_to_bop_dictionary(data):
    data = data.sort_values(by='MAPINFO')
    dic = {}
    for index in data.index:
        dic.update({data.at[index, 'ID_REF']: data.at[index, 'BOP']})
    return dic


def get_bop_to_cpg_dictionary(data):
    data = data.sort_values(by='BOP')
    dic = {}
    for index in data.index:
        bop = data.at[index, 'BOP']
        if bop not in dic:
            dic.update({bop: [data.at[index, 'ID_REF']]})
        else:
            dic[bop].append(data.at[index, 'ID_REF'])
    return dic


def main():
    try:
        annotations = open('..\\methylation_data\\annotations.dat', 'rb')
    except FileNotFoundError:
        print('There is no pickle file. Run create_table script.')
        return 2

    data = pickle.load(annotations)
    data = correct_data(data)
    print(data)

    cpg_to_bop = get_cpg_to_bop_dictionary(data)
    print(len(cpg_to_bop))
    i = 0
    for key, value in cpg_to_bop.items():
        print(key + ' - ' + str(value))
        i += 1
        if i > 9:
            break

    bop_to_cpg = get_bop_to_cpg_dictionary(data)
    print(len(bop_to_cpg))
    i = 0
    for key, value in bop_to_cpg.items():
        print(key + ' - ' + str(value))
        i += 1
        if i > 9:
            break

    result_file = open('result_bop_table.txt', 'w', encoding='utf-8')
    try:
        data_file = open('..\\methylation_data\\average_beta.txt', 'r', encoding='utf-8')
    except FileNotFoundError:
        print('There is no data file')
        return 2

    names = data_file.readline().split()
    width = len(names) - 1
    print('width =', width)
    table_dictionary = {}
    for key, value in bop_to_cpg.items():
        table_dictionary.update({key: np.zeros(width, dtype=np.float)})

    key_errors = 0
    for line in data_file:
        line = line.split()
        cpg_name = line.pop(0)
        line = np.array(line, dtype=np.float)
        try:
            bop = cpg_to_bop[cpg_name]
            table_dictionary[bop] += line
        except KeyError:
            key_errors += 1
        except ValueError:
            print(bop)
    print('KeyError count:', key_errors)

    for key, value in table_dictionary.items():
        value /= len(bop_to_cpg[key])
        result_file.write(key + '\t')
        for el in value:
            result_file.write(' ' + str(el))
        result_file.write('\n')

    result_file.close()
    data_file.close()
    annotations.close()


main()
