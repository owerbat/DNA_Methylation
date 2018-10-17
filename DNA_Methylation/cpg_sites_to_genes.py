import pandas
import pickle
import numpy as np


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
    data.UCSC_REFGENE_NAME = data.UCSC_REFGENE_NAME.apply(lambda string: remove_repeating_elements(string.split(';')))
    data.UCSC_REFGENE_ACCESSION = data.UCSC_REFGENE_ACCESSION.apply(lambda string:
                                                                    remove_repeating_elements(string.split(';')))
    data.UCSC_REFGENE_GROUP = data.UCSC_REFGENE_GROUP.apply(lambda string: remove_repeating_elements(string.split(';')))

    return data


def get_cpg_to_gene_dictionary(data):
    dic = {}
    for index in data.index:
        dic.update({data.at[index, 'ID_REF']: data.at[index, 'UCSC_REFGENE_NAME']})
    return dic


def main():
    try:
        file = open('annotations.dat', 'rb')
    except FileNotFoundError:
        print('There is no pickle file. Run create_table script.')
        return 2

    data = pickle.load(file)
    data = correct_data(data)
    print(data)
    print(list(data.columns))
    cpg_to_gene = get_cpg_to_gene_dictionary(data)
    print(len(cpg_to_gene))


main()
