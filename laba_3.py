import numpy as np
from laba_2 import main_phase


def first_phase(m, n, matrix_A, vector_b):
    # выполняем условие с вектором b
    for i in range(m):
        if vector_b[i] < 0:
            vector_b[i] *= (-1)
            matrix_A[i, :] *= (-1)
    # Объединяем А (входная матрица) и E (единичная матрица)
    matrix_AE = np.concatenate((matrix_A, np.eye(m)), axis=1)
    # небазисные индексы. Jn
    Jn = np.array(range(1, n + 1))
    # Базисные индексы. Jb
    vector_Jb_initial = np.array(range(n + 1, n + m + 1))
    # Родные переменные
    nativ_x = np.zeros(n)
    # Искусственные переменные
    not_nativ_x = np.copy(vector_b)
    # Новый комплект переменных x (c родными(данными) и искуственными переменными)
    new_vector_x = np.append(nativ_x, not_nativ_x)
    # Генерируем Вектор c
    vector_c = np.append(np.zeros(n), [-1] * m)
    # Оптимальный план для задачи и Jb
    # Решаем вспомогательную задачу основной фазой симплекс метода
    optimal_x, vector_Jb = main_phase(m, n, matrix_AE, None, vector_c, new_vector_x, vector_Jb_initial)
    optimal_x, vector_Jb = [round(i, 5) for i in optimal_x], list(vector_Jb)

    # если хотя бы одна искуственная переменная не равна 0, то задача несовместна и выходим из проги
    if not all(item == 0 for item in optimal_x[n:]):
        # Искуственные переменные не равны нулю. Задача несовместна.
        print("No solution")
        return None, None
    # Искуственные переменные равны нулю. Задача совместна. Продолжаем.
    # Находим индексы (jk), которые больше чем n и перебираем их
    # k = jk - n - 1
    for jk in [vector_Jb[i] for i in range(len(vector_Jb)) if vector_Jb[i] > n]:
        check = True
        for j in range(1, n + 1):
            if j in vector_Jb:
                check = True
                break
            # получаем базисную матрицу B
            matrix_B = np.zeros((m, len(vector_Jb)))
            for item in range(len(vector_Jb)):
                matrix_B[:, item] = matrix_A[:, vector_Jb[item] - 1]
            vector_l = np.matmul(np.linalg.inv(matrix_B), matrix_A[:, j - 1])
            if vector_l[jk - n - 1] != 0:
                vector_Jb[jk - n - 1] = j
                check = False
                continue
        if check == True:
            matrix_A = np.delete(matrix_A, jk - n - 1)
            vector_Jb.remove(vector_Jb[jk - n - 1])

    return optimal_x[:n], vector_Jb


if __name__ == '__main__':
    # variant 14 !!! (это он хочет из файла) 1.2.14
    m, n = 3, 4
    matrix_A = np.array([[-3, 7, 1, 0],
                  [7, 5, 0, 1],
                  [20, -4, -2, 2]])
    vector_b = np.array([14, 42, 56])

    # test 1 +
    '''
    m, n = 2, 2
    matrix_A = np.array([[1, 2],
                  [2, 1]])
    vector_b = np.array([3, 3])
    '''
    # test 2 +
    '''
    m, n = 1, 2
    matrix_A = np.array([[1, -1]])
    vector_b = np.array([0])
    '''
    # test 3
    '''
    m, n = 1, 2
    matrix_A = np.array([[1, 1]])
    vector_b = np.array([-1])
    '''

    x, j_b = first_phase(m, n, matrix_A, vector_b)
    print('\nResult:')
    # Базисный допустимый план X:
    print(x)
    # Множество базисных индексов, которое соответствует Х:
    print(j_b)

