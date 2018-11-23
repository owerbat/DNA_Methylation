from matplotlib import pyplot as plt
import numpy as np


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


def main():
    get_age_histogram('..\\methylation_data\\attributes GSE40279.txt')
    get_age_histogram('..\\methylation_data\\attributes GSE87571.txt')


main()
