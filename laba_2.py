import numpy as np
import math
from main import sherman_morrison_formula

# Реализация алгоритма.
# 1. Находим базисную матрицу Ab. Для этого составляем матрицу из нулевых элементов размерностью m и Jb.
#    После заполняем её столбцами, которые находятся под номерами базисных индексов.
# 2. Находим обратную базисную матрицу B.
# 3. Находим базисный вектор vector_cb.
# 4. Находим вектор потенциалов vector_u. Для этого умножаем базисный вектор vector_cb на обратную базисную матрицу B.
# 5. Находим вектор оценок vector_delta. Для этого умножаем вектор потенциалов vector_u на матрицу А и отнимаем вектор с
# 6. Проверяем что все компоненты вектора оценок vector_delta положетельные.
#    6.1. Если данное условие выполнено, то выводим сообщение 'Bounded' и вектор х. Выходим из программы.
#    6.2. Если данное условие не выполнено, и хотя бы один компонент вектора оценок vector_delta отрицательный,
#    то переходим к следующему пункту.
# 7. Находим минимальное значение min_num в векторе оценок vector_delta. Также находим индекс J0 этого минимального
#    значения.
# 8. Находим вектор z. Для этого умножаем обратную базисную матрицу B на столбец под номером J0 матрицы А.
# 9. Находим ϴ1, ..., ϴn. Для этого будем использовать вектор z.
#    9.1. Если i-тая компонента вектора z больше нуля то, вычисляем ϴi = vector_x[Jbi - 1] / vector_zi.
#    9.2. Если i-тая компонента вектора z меньше нуля то, ϴi равно бесконечности.
# 10. Находим минимальное значение min_teta_num из всех ϴ1, ..., ϴn..Также найдем значение индекса min_J этого элемента.
#    10.1. Если минимальное значение min_teta_num равно бесконечности, то оптимального плана нет. Выводим сообщение
#    "Unbounded" и выходим из программы.
#    10.2. Если минимальное значение min_teta_num не равно бесконечности, то переходим к следующему пункту.
# 11. Заменяем в исходном векторе Jb элемент, который стоит под индексом min_J, на индекс J0 минимального значения
#     в векторе оценок vector_delta (см. пункт 7).
# 12. Находим новый вектор x.
#    12.1. Вычисляем значения нового вектора х по формуле:
#    new_vector_xJi= vector_xJi - min_teta_num * vector_zi
#   12.2. Значение нового вектора х с индексом J0 равно min_teta_num.
# 13. Находим новое значение обратной базисной матрицы B. Повторяем алгоритм с пункта 3.


def positive_values(delta):
    val = True
    for item in delta:
        if item < 0:
            val = False
    return val


def min_and_index(value):
    item = min(np.array(value))
    for i in range(len(value)):
        if value[i] == item:
            item_index = i
            return item, item_index


def get_tetes(vector_z, vector_Jb, vector_x):
    len_z = len(vector_z)
    tetes = np.zeros(len_z)
    for i in range(len_z):
        if vector_z[i] > 0:
            tetes[i] = vector_x[vector_Jb[i] - 1] / vector_z[i]
        else:
            tetes[i] = math.inf
    return tetes


def get_vector_x_new(vector_x, vector_z, m, vector_Jb, min_teta, min_index_delta):
    vector_x_new = np.zeros(len(vector_x))
    for i in range(m):
        index = vector_Jb[i] - 1
        vector_x_new[index] = vector_x[index] - min_teta * vector_z[i]
    vector_x_new[min_index_delta] = min_teta
    return vector_x_new


def main_phase(m, n, matrix_A, vector_b, vector_c, vector_x, vector_Jb):
    # пунк 1
    matrix_Ab = np.zeros((m, len(vector_Jb)))
    for item in range(len(vector_Jb)):
        matrix_Ab[:, item] = matrix_A[:, vector_Jb[item] - 1]
    # пунк 2
    matrix_B = np.linalg.inv(matrix_Ab)
    while True:
        # пунк 3
        vector_cb = np.array([vector_c[item - 1] for item in vector_Jb])
        # пункт 4
        vector_u = np.dot(vector_cb, matrix_B)
        # пункт 5
        delta = np.dot(vector_u, matrix_A) - vector_c
        # пункт 6
        if positive_values(delta):
            print('Bounded')
            return vector_x, vector_Jb
        # пункт 7    J0 = min_index_delta
        min_delta, min_index_delta = min_and_index(delta)
        # пункт 8
        vector_z = np.dot(matrix_B, matrix_A[:, min_index_delta])
        # пункт 9
        tetes = get_tetes(vector_z, vector_Jb, vector_x)
        # пункт 10  (min_J = min_index_teta)
        min_teta, min_index_teta = min_and_index(tetes)
        if min_teta == math.inf:
            print("Unbounded")
            return None, None
        # пункт 11
        vector_Jb[min_index_teta] = min_index_delta + 1
        # пункт 12
        vector_x_update = get_vector_x_new(vector_x, vector_z, m, vector_Jb, min_teta, min_index_delta)
        vector_x = vector_x_update.copy()
        # пункт 13
        # calculations(B, A, vector_x, i, n):
        matrix_B = sherman_morrison_formula(m, min_index_teta + 1, matrix_Ab, matrix_B, matrix_A[:, min_index_delta])
        if matrix_B is None:
            print("Unbounded")
            return None, None


if __name__ == '__main__':
    # test 1 +
    '''
    m, n = 2, 2
    matrix_A = np.array([[1, -1],
                         [2, 1]])
    vector_b = np.array([3, 3])
    vector_c = np.array([1, 1])
    vector_x = np.array([1, 1])
    vector_Jb = [2, 1]
    '''
    # test 2 +
    '''
    m, n = 1, 2
    matrix_A = np.array([[1, -1]])
    vector_b = np.array([0])
    vector_c = np.array([1, 2])
    vector_x = np.array([0, 0])
    vector_Jb = [1]
    '''

    m, n = 4, 6
    matrix_A = np.array([[5, 2, 1, 0, 0, 0],
                         [2, 7, 0, 1, 0, 0],
                         [1, 0, 0, 0, 1, 0],
                         [0, 1, 0, 0, 0, 1]])
    vector_b = np.array([45, 49, 8, 6])
    vector_c = np.array([25, 32, 1, 0, 0, 0])
    vector_x = np.array([0, 0, 45, 49, 8, 6])
    vector_Jb = [3, 4, 5, 6]

    result_vector_x, result_vector_Jb = main_phase(m, n, matrix_A, vector_b, vector_c, vector_x, vector_Jb)
    print(result_vector_x)
