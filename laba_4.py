import numpy as np
import math


def sherman_morrison_formula(n, i, A, B, vector_x):
    vector_l = np.dot(B, vector_x)
    if vector_l[i] != 0:
        vector_l_wave = vector_l.copy()
        vector_l_wave[i] = -1
        vector_l_cap = vector_l_wave.dot(((-1) / vector_l[i]))
        matrix_E = np.eye(n)
        matrix_E[:, i] = vector_l_cap
        Q = np.dot(matrix_E, B)
        return Q
    else:
        return None


def min_and_index(value):
    item = min(np.array(value))
    for i in range(len(value)):
        if value[i] == item:
            item_index = i
            return item, item_index


def positive_values(value):
    res = True
    for item in value:
        if item < 0:
            res = False
    return res


def dual_simplex_method(c, b, a_matrix, j_vector, m, n):
    # базисня матрица
    Ab = a_matrix[:, j_vector]
    # обратная базисная матрица
    B = np.linalg.inv(Ab)
    # вектор c_b (берем базисные индексы из начального вектора c)
    cb = np.array([c[j] for j in j_vector])
    # начальный баз. двойственный план
    y = np.dot(cb, np.linalg.inv(Ab))
    # дальше итерации
    while True:
        # находим вектор Jn (значения, которые не входят в вектор Jb)
        Jn = [j for j in range(n) if j not in j_vector]
        # вычисляем базисные индексы каппа (перед началом создали вектор x_0 т.к. небазисные индексы вектора каппа = 0)
        vector_kappa = np.dot(B, b)
        vector_kappa_0 = [0 for _ in range(n)]
        # если все значения вектора каппа положительные, то решение найдено и выводим результат
        if positive_values(vector_kappa):
            # объединяем базисные и небазисные индесы вектора kappa
            for j, kappa_item in zip(j_vector, vector_kappa):
                vector_kappa_0[j] = kappa_item
            return vector_kappa_0
        # находим минимальное значение вектора kappa (kappa_min)
        # находим индекс минимальной компоненты вектора kappa (kappa_min_index)
        kappa_min, kappa_min_index = min_and_index(vector_kappa)
        # берем строку B под индексом минимального компонента
        delta_y = B[kappa_min_index]
        # находим вектор mu
        mu = delta_y.dot(a_matrix)
        # если все mu положительные, то прямая задача несовместна (нет допустимого плана)
        if positive_values(mu):
            print("No solution")
            return None
        sigma = []
        # для всех Jn(j) для которых mu(j) > 0 находим вектор sigma
        for i in Jn:
            if mu[i] < 0:
                sigma.append((c[i] - np.dot(a_matrix[:, i], y)) / mu[i])
            else:
                sigma.append(math.inf)
        # находим минимальное значения из всех сигм (sigma_min)
        # находи индекс минимального значения (запоминаем его) (sigma_min_index)
        sigma_min, sigma_min_index = min_and_index(sigma)
        if sigma_min == math.inf:
            print("Задача не имеет решения, т.к. пусто множество ее допустимых планов.")
            return None
        # заменяем Jb под индексом минимальной компоненты вектора kappa (kappa_min_index) на индекс
        # минимального значения вектора сигмы (sigma_min_index)
        j_vector[kappa_min_index] = sigma_min_index
        # обновим y
        y = y + np.dot(sigma_min, delta_y)
        # обновим Ab
        Ab = a_matrix[:, j_vector]
        # обновим B (используем 1-ую лабу)
        B = sherman_morrison_formula(m, kappa_min_index, Ab, B, A[:, sigma_min_index])


if __name__ == "__main__":
    m, n = 2, 4
    A = np.array([[1, -5, 1, 0],
                  [0, -4, 1, 1]])
    c = np.array([-2, 0, 1, 0])
    b = np.array([-10, -12])
    J = np.array([2, 3])
    # Ответ -> [5, 3, 0, 0] (вариант 14)

    '''
    m, n = 2, 4
    A = np.array([[2, -3, 1, 0],
                  [-1, -1, 1, 1]])
    c = np.array([1, -4, 1, 0])
    b = np.array([-3, -6])
    J = np.array([2, 3])
    '''
    # Ответ -> [3, 3, 0, 0] вариант 12
    x = dual_simplex_method(c, b, A, J, m, n)
    print(x)