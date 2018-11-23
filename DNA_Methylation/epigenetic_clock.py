from matplotlib import pyplot as plt
import numpy as np


def get_age_histogram():
    file = open('..\\methylation_data\\attributes GSE40279.txt', 'r', encoding='utf-8')
    ages = []
    male_ages = []
    female_ages = []
    file.readline()
    for line in file:
        age = int(line.split()[2])
        ages.append(age)

        sex = line.split()[3]
        if sex == 'M':
            male_ages.append(age)
        else:
            female_ages.append(age)

    plt.xlabel('Age')
    plt.ylabel('Count')
    plt.hist(ages, bins=np.arange(min(ages)+0.5, max(ages), 1), density=True, histtype='step', color='black')
    plt.hist(male_ages, bins=np.arange(min(ages) + 0.5, max(ages), 1), density=True, histtype='step', color='blue')
    plt.hist(female_ages, bins=np.arange(min(ages) + 0.5, max(ages), 1), density=True, histtype='step', color='pink')
    plt.show()


def main():
    get_age_histogram()


main()
