import math
import argparse
from functools import reduce


def st(a, b, n):
    m = max(n) + 1
    return [x / m * b + a for x in n]


def tr(a, b, n):
    m = max(n) + 1
    n = [x / m for x in n]
    return [a + b * (n[i] + n[i + 1] - 1) for i in range(len(n) - 1)]


def ex(a, b, n):
    m = max(n) + 1
    n = [x / m for x in n]
    return [-b * math.log(u) + a for u in n]


def nr(a, b, n):
    result = []
    m = max(n) + 1
    n = [x / m for x in n]
    for i in range(0, len(n) - 1, 2):
        z1 = a + b * math.sqrt(-2 * math.log(1 - n[i])) * math.cos(2 * math.pi * n[i + 1])
        z2 = a + b * math.sqrt(-2 * math.log(1 - n[i])) * math.sin(2 * math.pi * n[i + 1])
        result.append(z1)
        result.append(z2)
    return result


def gm(a, b, k, n):
    m = max(n) + 1
    n = [x / m for x in n]
    list = [n[i: i + k] for i in range(0, len(n), k)]
    if len(list[-1]) != k:
        list.pop()
    return [a - b * math.log(reduce(lambda x, y: x * (1 - y), first, 1)) for first in list]


def ln(a, b, n):
    n = nr(0, 1, n)
    return [a + math.exp(b - z) for z in n]


def ls(a, b, n):
    m = max(n) + 1
    n = [x / m for x in n]
    return [a + b * math.log(u / (1 - u)) for u in n]


def bi(a, num):
    n = len(num)
    m = max(num) + 1
    num = [x / m for x in num]
    result = []
    for i in range(1, n):
        temp = 0
        y = num[i]
        while temp < num[i]:
            k = 0
            while y > k:
                temp += math.comb(n, k) * pow(a, k) * pow(1 - a, n - k)
                k += 1
            y += 1
        result.append(y * 10000)
    return result


def write_to_file(data, filepath):
    with open(filepath, "w", encoding="UTF-8") as f:
        f.write(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="rnc.exe",
        description="Генерация последовательности псевдослучайных чисел по выбранному распределению")
    help1 = """Указывает тип распределения: 
                st – стандартное равномерное с заданным интервалом,
                tr – треугольное распределение,
                ex – общее экспоненциальное распределение, 
                nr – нормальное распределение,
	            gm – гамма распределение,
	            ln – логнормальное распределение, 
	            ls – логистическое распределение,
	            bi – биномиальное распределение"""
    help2 = "Указывает путь до файла из которого берется входная последовательность чисел"
    help3 = "1-й параметр, необходимый, для генерации ПСЧ заданного распределения"
    help4 = "2-й параметр, необходимый, для генерации ПСЧ заданного распределения (в случае bi указывает на количество разрядов в генерируемом чиле до 10^6)"
    help5 = "3-й параметр, необходимый, для генерации ПСЧ гамма-распределением."
    parser.add_argument("-d", help=help1, required=True, choices=["st", "tr", "ex", "nr", "gm", "ln", "ls", "bi"],
                        nargs=1)
    parser.add_argument("-f", nargs=1, default=["rnd.txt"], help=help2)
    parser.add_argument("-p1", nargs=1, type=int, required=True, help=help3)
    parser.add_argument("-p2", nargs=1, type=int, help=help4, default=[None])
    parser.add_argument("-p3", nargs=1, type=int, help=help5, default=[None])
    args = parser.parse_args()


    def check_params(count_params, parametrs):
        if len(parametrs) != count_params:
            raise Exception("Передано неверное количество аргументов")
        for par in parametrs:
            if not par.isdigit():
                raise Exception("Переданы неподходящие параметры")
        return True


    type = args.d[0]
    path_file = args.f[0]
    p1 = args.p1[0]
    p2 = args.p2[0]
    p3 = args.p3[0]
    try:
        with open(path_file, "r") as f:
            line = f.readline()
            list = list(map(int, line.split(",")))
        match type:
            case 'st':
                a = p1
                b = p2
                file = "distr-st.dat"
                numbers = st(a, b, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'tr':
                a = p1
                b = p2
                file = "distr-tr.dat"
                numbers = tr(a, b, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'ex':
                a = p1
                b = p2
                file = "distr-ex.dat"
                numbers = ex(a, b, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'nr':
                a = p1
                b = p2
                file = "distr-nr.dat"
                numbers = nr(a, b, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'gm':
                a = p1
                b = p2
                k = p3
                file = "distr-gm.dat"
                numbers = gm(a, b, k, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'ln':
                a = p1
                b = p2
                file = "distr-ln.dat"
                numbers = ln(a, b, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'ls':
                a = p1
                b = p2
                file = "distr-ls.dat"
                numbers = ls(a, b, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
            case 'bi':
                a = p1
                file = "distr-bi.dat"
                numbers = bi(a, list)
                data = ",".join(map(str, numbers))
                write_to_file(data, filepath=file)
    except Exception as err:
        print("В процессе генерации произошла ошибка!")
        print(" • " + str(err))
