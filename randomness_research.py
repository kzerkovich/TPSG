import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from math import factorial, floor
import argparse


def mean(seq):
    return np.mean(seq)


def relative_error_mean(seq):
    mean_val = mean(seq)
    std_dev = np.std(seq, ddof=1)
    relative_error = std_dev / np.sqrt(len(seq)) / mean_val
    return relative_error


def std_dev(seq):
    return np.std(seq, ddof=1)


def relative_error_std_dev(seq):
    std_dev_val = std_dev(seq)
    n = len(seq)
    relative_error = std_dev_val / np.sqrt(2 * (n - 1))
    return relative_error


def normalize(seq):
    min_val = np.min(seq)
    max_val = np.max(seq)
    return (seq - min_val) / (max_val - min_val + 1)


def plot_statistics(seq, step=500):
    normalized_seq = normalize(abs(seq)) if np.max(seq) > 1 else seq
    sample_sizes = list(range(step, len(normalized_seq) + 1, step))
    means = []
    std_devs = []

    for size in sample_sizes:
        sample = normalized_seq[:size]
        means.append(mean(sample))
        std_devs.append(std_dev(sample))

    plt.figure(figsize=(14, 7))

    plt.subplot(1, 2, 1)
    plt.plot(sample_sizes, means, label='Математическое ожидание', marker='o')
    plt.xlabel('Длина выборки')
    plt.ylabel('Математическое ожидание')
    plt.title('Зависимость математического ожидания от выборки')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(sample_sizes, std_devs, label='Среднеквадратичное отклонение', marker='o')
    plt.xlabel('Длина выборки')
    plt.ylabel('Среднеквадратичное отклонение')
    plt.title('Зависимость среднеквадратичного отклонения от выборки')
    plt.legend()

    plt.tight_layout()
    plt.show()


def count(seq):
    lst = {}
    for x in seq:
        if x in lst:
            lst[x] += 1
        else:
            lst[x] = 1
    return lst


def xi_square(seq, num_list=None, k=None, e=None, a=0.05):
    if num_list is None:
        num_list = list(count(seq).values())
    if k is None:
        k = len(num_list)
    if e is None:
        e = [len(seq) / k] * k
    xi = sum(map(lambda x: ((x[0] - x[1]) ** 2) / x[1], zip(num_list, e)))
    squa = chi2.ppf(1 - a, k - 1)
    return xi <= squa


def series(seq):
    d = 64
    a = 0.05
    res = [0] * (d ** 2)
    for j in range(len(seq) // 2):
        q = int(seq[2 * j] * d)
        r = int(seq[2 * j + 1] * d)
        res[q * d + r] += 1
    e = d ** 2 * [len(seq) / (2 * d ** 2)]
    return xi_square(seq=seq, num_list=res, k=d * d, e=e, a=a)


def interval(seq):
    a = 0.05
    t = 10
    f = 0.5
    b = 1.0
    j = -1
    s = 0
    c = t * [0]
    n = len(seq)
    intervals = n / 10
    while s < intervals and j < n:
        r, j = 0, j + 1
        while f <= seq[j] <= b and j < n:
            r, j = r + 1, j + 1
        c[min(r, t) - 1] += 1
        s += 1
    sub = b - f
    e = [intervals * sub * (1 - sub) ** z for z in range(t)]
    return xi_square(seq=seq, num_list=c, e=e, k=t + 1, a=a)


def stirling(r, k):
    if r <= 0 or r != 0 and r == k:
        return 1
    elif k <= 0 or r < k:
        return 0
    elif r == 0 and k == 0:
        return - 1
    else:
        return k * (stirling(r - 1, k)) + stirling(r - 1, k - 1)


def partition(seq):
    d = 64
    k = 5
    n = len(seq)
    num_list = [0] * k
    for i in range(n // k):
        unique = list(set([int(u * d) for u in seq[i * k: i * k + k]]))
        num_list[len(unique) - 1] += 1
    e = [0] * k
    for r in range(1, k + 1):
        p = 1.0
        for i in range(r):
            p *= d - i
        e[r - 1] = (n / k) * (p / d ** k) * stirling(k, r)
    return xi_square(seq=seq, num_list=num_list, k=k, e=e)


def permutation(seq):
    t = 10
    a = 0.05
    n = len(seq)
    k = factorial(t)
    tmp = {}
    for i in range(0, n, t):
        group = tuple(sorted(seq[i:i + t]))
        tmp[group] = tmp.get(group, 0) + 1
    num_list = list(tmp.values())
    e = [n / k] * len(num_list)
    return xi_square(seq=seq, num_list=num_list, k=k, e=e, a=a)


def monoton(seq):
    a = 0.05
    matrix = [[4529.4, 9044.9, 13568.0, 18091.0, 22615.0, 27892.0],
              [9044.9, 18097.0, 27139.0, 36187.0, 45234.0, 55789.0],
              [13568.0, 27139.0, 40721.0, 54281.0, 67852.0, 83685.0],
              [18091.0, 36187.0, 54281.0, 72414.0, 90470.0, 111580.0],
              [22615.0, 45234.0, 67852.0, 90470.0, 113262.0, 139476.0],
              [27892.0, 55789.0, 83685.0, 111580.0, 139476.0, 172860.0]]
    vector = [1.0 / 6.0, 5.0 / 24.0, 11.0 / 120.0, 19.0 / 720.0, 29.0 / 5040.0, 1.0 / 840.0]
    n = len(seq)
    tmplst = []
    i = 0
    while i < n:
        s = 1
        while i + s < n and seq[i + s - 1] <= seq[i + s]:
            s += 1
        tmplst.append(s)
        i += s
    m = 0
    group = {}
    for length in tmplst:
        m = max(m, length)
        group[length] = group.get(length, 0) + 1
    e = []
    tmp = 0
    for c in tmplst:
        m = 1 / 6
        min_value = min(c, 6)
        for i in range(min_value):
            for j in range(min_value):
                m += (seq[i + tmp] - n * vector[i]) * (seq[j + tmp] - n * vector[j]) * matrix[i][j]
        tmp += c
        e.append(m)
    return xi_square(seq=seq, e=e, a=a)


def check_conflicts(seq):
    n = len(seq)
    m = pow(2, 10)
    s = n / m
    p0 = 1 - n / m + factorial(n) / (2 * factorial(n - 2) * pow(m, 2))
    conf = n / m - 1 + p0
    return conf < s + 10 or conf > s - 10


def get_data_from_file(file_path):
    with open(file_path, "r") as file:
        data_line = file.readline()
        num_lst = list(map(int, data_line.strip().split(",")))
    return normalize(num_lst)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="check.py",
        description="Проверка последовательности чисел по заданному критерию"
    )
    help1 = "Имя файла с входной последовательностью"
    help2 = """ Критерии для проверки заданной последовательности на равномерное распределение:
                  •  chi - Критерий хи-квадрат,
                  •  series - Критерий серий,
                  •  intervals - Критерий интервалов,
                  •  partition - Критерий разбиений,
                  •  permutation - Критерий перестановок,
                  •  monoton - Критерий монотонности,
                  •  conflict - Критерий конфликтов
         """
    types_crit = ["chi", "series", "intervals", "partition", "permutation", "monoton", "conflict"]
    parser.add_argument("-type", nargs=1, choices=types_crit, required=False, default=[None], help=help2)
    parser.add_argument("-file", nargs=1, required=False, default=[None], help=help1)
    args = parser.parse_args()
    type_criteries = args.type[0]
    file_path = args.file[0]
    try:
        seq = get_data_from_file(file_path)
        num1 = mean(seq)
        num2 = std_dev(seq)
        num3 = relative_error_mean(seq)
        num4 = relative_error_std_dev(seq)
        print(f"Мат. ожидание: {num1}")
        print(f"Cр-кв. отклонение: {num2}")
        print(f"Относительное мат. отклонение: {num3}")
        print(f"Относительное ср-кв. отклонение: {num4}")
        if type_criteries == None:
            print("Критерий хи-квадрат:")
            if (xi_square(seq)):
                print("• Удовлетворяет критерию хи-квадрат")
            else:
                print("• Не удовлетворяет критерию хи-квадрат")
            print("Критерий серий:")
            if series(seq):
                print("• Удовлетворяет критерию серий")
            else:
                print("• Не удовлетворяет критерию серий")
            print("Критерий интервалов:")
            if interval(seq):
                print("• Удовлетворяет критерию интервалов")
            else:
                print("• Не удовлетворяет критерию интервалов")
            print("Критерий разбиений:")
            if partition(seq):
                print("• Удовлетворяет критерию разбиений")
            else:
                print("• Не удовлетворяет критерию разбиений")
            print("Критерий перестановок:")
            if permutation(seq):
                print("• Удовлетворяет критерию перестановок")
            else:
                print("• Не удовлетворяет критерию перестановок")
            print("Критерий монотонности:")
            if monoton(seq):
                print("• Удовлетворяет критерию монотонности")
            else:
                print("• Не удовлетворяет критерию монотонности")
            print("Критерий конфликтов:")
            if check_conflicts(seq):
                print("• Удовлетворяет критерию конфликтов")
            else:
                print("• Не удовлетворяет критерию конфликтов")
        else:
            match type_criteries:
                case 'chi':
                    print("Критерий хи-квадрат:")
                    if (xi_square(seq)):
                        print("• Удовлетворяет критерию хи-квадрат")
                    else:
                        print("• Не удовлетворяет критерию хи-квадрат")
                case 'series':
                    print("Критерий серий:")
                    if series(seq):
                        print("• Удовлетворяет критерию серий")
                    else:
                        print("• Не удовлетворяет критерию серий")
                case 'intervals':
                    print("Критерий интервалов:")
                    if interval(seq):
                        print("• Удовлетворяет критерию интервалов")
                    else:
                        print("• Не удовлетворяет критерию интервалов")
                case 'partition':
                    print("Критерий разбиений:")
                    if partition(seq):
                        print("• Удовлетворяет критерию разбиений")
                    else:
                        print("• Не удовлетворяет критерию разбиений")
                case 'permutation':
                    print("Критерий перестановок:")
                    if permutation(seq):
                        print("• Удовлетворяет критерию перестановок")
                    else:
                        print("• Не удовлетворяет критерию перестановок")
                case 'monoton':
                    print("Критерий монотонности:")
                    if monoton(seq):
                        print("• Удовлетворяет критерию монотонности")
                    else:
                        print("• Не удовлетворяет критерию монотонности")
                case 'conflict':
                    print("Критерий конфликтов:")
                    if check_conflicts(seq):
                        print("• Удовлетворяет критерию конфликтов")
                    else:
                        print("• Не удовлетворяет критерию конфликтов")
    except Exception as err:
        print("В процессе проверки произошла ошибка!")
        print(" • " + str(err))
    plot_statistics(seq)
