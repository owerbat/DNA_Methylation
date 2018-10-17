import pandas
import pickle


def get_data_frame(file_name):
    try:
        file = open(file_name)
    except FileNotFoundError:
        print('The file is not founded')

    names = file.readline().split()
    data_dict = {}
    i = 0
    while True:
        current = file.readline().split('\t')
        if current == ['']:
            break
        data_dict.update({i: current})
        i += 1
    data = pandas.DataFrame(data_dict, index=names)
    data = data.transpose()
    return data


def main():
    data = get_data_frame('D:\\DNA Methylation\\annotations.txt')
    file = open('annotations.dat', 'wb')
    pickle.dump(data, file)
    file.close()


main()
