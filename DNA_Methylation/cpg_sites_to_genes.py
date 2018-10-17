import pandas
import pickle


def correct_data(data):
    names = list(data.columns)
    names.pop(3)
    names.pop(3)
    data = data[names]
    data = data[data.RELATION_TO_UCSC_CPG_ISLAND != 'N_Shelf']
    data = data[data.RELATION_TO_UCSC_CPG_ISLAND != 'S_Shelf']
    data = data[data.CHR != 'X']
    data = data[data.CHR != 'Y']
    return data


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


main()
