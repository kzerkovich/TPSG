import argparse


def progress_bar(now, all):
    if ((now + 1) / all * 100) % 20 == 0:
        print(f"\nГенерация: {(now + 1) / all * 100}%")


def lc(m, a, c, x0, n=10000):
    res = [x0]
    for i in range(n):
        res.append((a * res[i] + c) % m)
        progress_bar(i, n)

    return res


def add(m, li, ui, arr, n=10000):
    res = []
    tmp = arr.copy()
    for i in range(n):
        new_elem = (tmp[li - 1] + tmp[ui - 1]) % m
        res.append(new_elem)
        tmp.append(new_elem)
        tmp.pop(0)
        progress_bar(i, n)
    return res


def p5(p, q1, q2, q3, w, x0, n=10000):
    res = []
    x0 = int(x0, base=2)
    for i in range(n):
        new_elem = 1
        for _ in range(w - 1):
            new_bit1 = (x0 >> p - q1) & 1
            new_bit2 = (x0 >> p - q2) & 1
            new_bit3 = (x0 >> p - q3) & 1
            last_bit = p & 1
            xor = new_bit1 ^ new_bit2 ^ new_bit3 ^ last_bit
            new_elem = (new_elem << 1) | xor
            x0 = (x0 >> 1) | (xor << p - 1)
        progress_bar(i, n)
        res.append(new_elem)
    return res


def lfsr(vec_of_cf, register, n=10000):
    res = []
    len_reg = len(register)
    vec_of_cf = int(vec_of_cf, 2)
    register = int(register, 2)
    for i in range(n):
        new_bit = (register ^ ((register ^ vec_of_cf) >> 1)) & 1
        register = (register >> 1) | (new_bit << (len_reg - 1))
        res.append(register)
        progress_bar(i, n)
    return res


def nfsr(R1, R2, R3, w, x1, x2, x3, n=10000):
    res = []
    len_R1 = len(R1)
    len_R2 = len(R2)
    len_R3 = len(R3)
    for i in range(n):
        new_elem = 1
        for _ in range(w - 1):
            xor_R1 = (x1 ^ (x1 >> 1))
            xor_R2 = (x2 ^ (x2 >> 1))
            xor_R3 = (x3 ^ (x3 >> 1))
            out = ((xor_R1 ^ xor_R2) + (xor_R2 ^ xor_R3) + xor_R3) & 1
            x1 = (x1 >> 1) | (out << len_R1 - 1)
            x2 = (x2 >> 1) | (out << len_R2 - 1)
            x3 = (x3 >> 1) | (out << len_R3 - 1)
            new_elem = (new_elem << 1) | out
        res.append(new_elem)
        progress_bar(i, n)
    return res


def mt(x0, m, n=10000):
    res = []
    w = 32
    start = [0] * m
    start[0] = x0
    for i in range(1, m):
        start[i] = (start[i - 1] ^ (start[i - 1] >> 30) + i) & ((1 << w) - 1)
    for i in range(n):
        index = m
        u, s, t, l, w = 11, 7, 15, 18, 32
        b = 0x9D2C5680
        c = 0xEFC60000
        if index >= 624:
            p, r, q = 624, 31, 397
            a = 0x9908B0DF
            for j in range(p):
                rand_num = (start[j] >> r) + (start[(j + 1) % p] & ((1 << r) - 1))
                start[j] = start[(j + q) % p] ^ (rand_num >> 1)
                if rand_num % 2 != 0:
                    start[j] ^= a
            index = 0
        rand_num = start[index]
        rand_num ^= (rand_num >> u)
        rand_num ^= ((rand_num << s) & b)
        rand_num ^= ((rand_num << t) & c)
        rand_num ^= (rand_num >> l)
        index += 1
        rand_num = rand_num & ((1 << w) - 1)
        res.append(rand_num)
        progress_bar(i, n)
    return res


def rc4(vec, n=10000):
    vec_size = len(vec)
    S = list(range(4096))
    j = 0
    for i in range(4096):
        j = (j + S[i] + vec[i % vec_size]) % 4096
        S[i], S[j] = S[j], S[i]

    f = 0
    s = 0
    res = []
    for i in range(n):
        f = (f + 1) % 4096
        s = (s + S[f]) % 4096
        S[f], S[s] = S[s], S[f]
        new_elem = S[(S[f] + S[s]) % 4096]
        res.append(new_elem)
        progress_bar(i, n)
    return res


def rsa(n, e, w, x, num=10000):
    res = []
    for i in range(num):
        new_elem = 0
        for _ in range(w):
            x = pow(x, e, n)
            new_elem = (new_elem << 1) | (x & 1)
        res.append(new_elem)
        progress_bar(i, num)
    return res


def bbs(x, w, num=10000):
    p = 127
    q = 131
    n = p * q
    res = []
    for i in range(num):
        new_elem = 0
        for _ in range(w):
            x = pow(x, 2, n)
            new_elem = (new_elem << 1) | (x & 1)
        res.append(new_elem)
        progress_bar(i, num)
    return res


def write_to_file(data, filepath="rnd.dat"):
    with open(filepath, "w", encoding="UTF-8") as f:
        f.write(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="prng.exe",
        description="Генерация последовательности псевдослучайных чисел")
    help1 = """ - метод генерации:\n                        
            •   lc – линейный конгруэнтный метод;\n
            •   add – аддитивный метод;\n
            •   5p – пятипараметрический метод;\n
            •   lfsr – регистр сдвига с обратной связью (РСЛОС);\n
            •   nfsr – нелинейная комбинация РСЛОС;\n
            •   mt – вихрь Мерсенна;\n
            •   rc4 – RC4;\n
            •   rsa – ГПСЧ на основе RSA;\n
            •   bbs – алгоритм Блюма-Блюма-Шуба;\n
            """
    help2 = "- количество генерируемых чисел. По умолчанию - 10000."
    help3 = "- имя файла, в который будут записываться данные.\n По умолчанию - rnd.dat."
    help4 = "Перечисление параметров для выбранного генератора"
    parser.add_argument("-g", help=help1, required=True,
                        choices=["lc", "add", "5p", "lfsr", "nfsr", "mt", "rc4", "rsa", "bbs"], nargs=1)
    parser.add_argument("-i", nargs="*", help=help4, )
    parser.add_argument("-n", nargs=1, type=int, default=[10000], help=help2)
    parser.add_argument("-f", nargs=1, default=["rnd.dat"], help=help3, )
    args = parser.parse_args()


    def check_params(count_params, parametrs):
        if len(parametrs) != count_params:
            raise Exception("Передано неверное количество аргументов")
        for par in parametrs:
            if not par.isdigit():
                raise Exception("Переданы неподходящие параметры")
        return True


    type = args.g[0]
    parametrs = args.i
    path_file = args.f[0]
    count = args.n[0]
    try:
        match type:
            case 'lc':
                if check_params(4, parametrs):
                    m = int(parametrs[0])
                    a = int(parametrs[1])
                    c = int(parametrs[1])
                    x0 = int(parametrs[2])
                    result = lc(m, a, c, x0, n=count)
                data = ",".join(map(str, result))
                write_to_file(data, filepath=path_file)
            case 'add':
                if len(parametrs) < 4:
                    raise Exception("Неверное количество аргументов")
                for par in parametrs:
                    if not par.isdigit():
                        raise Exception("Неподходящие параметры")
                m = int(parametrs[0])
                li = int(parametrs[1])
                ui = int(parametrs[2])
                arr = list(map(int, parametrs[3:]))
                result = add(m, li, ui, arr, n=count)
                data = ",".join(map(str, result))
                write_to_file(data, filepath=path_file)
            case '5p':
                if check_params(6, parametrs):
                    p = int(parametrs[0])
                    q1 = int(parametrs[1])
                    q2 = int(parametrs[2])
                    q3 = int(parametrs[3])
                    w = int(parametrs[4])
                    x0 = parametrs[5]
                    result = p5(p, q1, q2, q3, w, x0, n=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
            case 'lfsr':
                if check_params(2, parametrs):
                    vec_of_cf = parametrs[0]
                    registr = parametrs[1]
                    result = lfsr(vec_of_cf, registr, n=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
            case 'nfsr':
                if check_params(7, parametrs):
                    R1 = parametrs[0]
                    R2 = parametrs[1]
                    R3 = parametrs[2]
                    w = int(parametrs[3])
                    x1 = int(parametrs[4])
                    x2 = int(parametrs[5])
                    x3 = int(parametrs[6])
                    result = nfsr(R1, R2, R3, w, x1, x2, x3, n=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
            case 'mt':
                if check_params(2, parametrs):
                    x0 = int(parametrs[1])
                    m = int(parametrs[0])
                    result = mt(x0, m, n=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
            case 'rc4':
                if check_params(256, parametrs):
                    vec = list(map(int, parametrs))
                    result = rc4(vec, n=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
            case 'rsa':
                if 4 == len(parametrs):
                    n = int(parametrs[0])
                    e = int(parametrs[1])
                    w = int(parametrs[2])
                    x = int(float(parametrs[3]))
                    result = rsa(n, e, w, x, num=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
            case 'bbs':
                if check_params(2, parametrs):
                    x = int(parametrs[0])
                    w = int(parametrs[1])
                    result = bbs(x, w, num=count)
                    data = ",".join(map(str, result))
                    write_to_file(data, filepath=path_file)
    except Exception as err:
        print("Произошла ошибка!")
        print(" • " + str(err))
