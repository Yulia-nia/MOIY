import numpy as np
import math


# условие баланса
def balance_condition(a, b, c):
    # находим сумму элементов a и b
    sum_a, sum_b = sum(a), sum(b)
    # сравниваем суммы, прибаляем разность сумм к элементам b и а
    # + добавляем к с
    if sum_a > sum_b:
        b.append(sum_a - sum_b)
        # добавляем к С столбец из 0
        c = np.append(c, np.zeros((len(a), 1)), axis=1)
    elif sum_a < sum_b:
        a.append(sum_b - sum_a)
        # добавляем к С строку из 0
        c = np.append(c, np.zeros((1, len(b))), axis=0)
    # елсли суммы равны, возвращаем неизмененные значения
    return a, b, c


# метод северо-западного угла
# Суть метода:
# Начинаем с первых индексов (0, 0). Для a хар-но i. Для b хар-но j.
# на каждой итерации определяем, сколько можно пренести из a в b. И записываем кол-во в X под индексами (i, j).
# если в a остаеться 0 то увеличиваем i на 1. если в b остаеться 0 то увеличиваем j на 1.
# (за один шаг можно увеличить либо i, либо j)
# и т.д.
# Можем сдвинуть i? Нет. Можем сдвинуть j? Нет. Звершаем работу.
def northwest_angle_method(a, b):
    m, n = len(a), len(b)
    # сначала план x -> нулевая матрица и вектор Jb -> пустой
    matrix_x = np.zeros((m, n))
    Jb = []
    i, j = 0, 0
    while True:
        Jb.append((i, j))
        max_supply = min(a[i], b[j])
        a[i] -= max_supply
        b[j] -= max_supply
        matrix_x[i][j] = max_supply
        if i == m - 1 and j == n - 1:
            # X - базисный план перевозок, Jb - множество базисных позиций.
            return matrix_x, Jb
        if a[i] == 0 and i != len(a) - 1:
            i += 1
        elif b[j] == 0 and j != len(b) - 1:
            j += 1


# находим U и V
def potentials_method(c, Jb):
    m, n = np.shape(c)
    a = np.zeros((m + n, m + n))
    b = np.zeros(m + n)
    # приравниваем u0 = 0 ( добавляем строчку вниз)
    count = 0
    a[-1][0] = 1
    b[-1] = 0
    # создаем матрицу для решения системы
    for i, j in Jb:
        a[count][i] = 1
        a[count][j + m] = 1
        b[count] = c[i][j]
        count += 1
    # решаем систему
    x = np.linalg.solve(a, b)
    u = x[:m]
    v = x[m:]
    return u, v


# проверяем условие оптимальности
# для любой небазисной позиции (i, j) выполняется условие -> Ui - Vj <= Cij
# елсли для какой-то позиции условие не выполняется (т.е. u[i] + v[j] > c[i][j]) --> находим индексы этого элемента
def check_optimality_condition(u, v, c, J_b):
    for i in range(np.shape(c)[0]):
        for j in range(np.shape(c)[1]):
            if (i, j) not in J_b:
                if u[i] + v[j] > c[i][j]:
                    # print('Условие оптимальности не выполняется в позиции ' + str((i, j)))
                    return i, j
    # отправляем что условие выполняется для всех (ведет к завершению программы)
    return None, None


# Проходимся по всем строкам матрицы X. Если в строке базисных позиций нет ли только одна, то удаляем её.
# Потом также по столбцам. Потом снова по строкам. И так пока мы ничего не можем удалить.
# Все оставшиеся базисные позиции - угловые.
def find_matrix_corner(m, n, Jb, I, J, x):
    Jb.append((I, J))
    no_lines_to_delete, no_rows_to_delete = False, False
    # временный Jb
    Jb_temp = Jb.copy()
    m_x, n_x = np.shape(x)
    while no_lines_to_delete is False or no_rows_to_delete is False:
        no_lines_to_delete = True
        no_rows_to_delete = True
        #m_x, n_x = np.shape(x)
        # строки
        # ищем в строке все базисные позиции
        for i in range(m_x):
            basis_indexes_count = 0
            for j in range(n_x):
                if (i, j) in Jb_temp:
                    basis_indexes_count += 1
            # удаляем
            if basis_indexes_count <= 1:
                for index in range(m):
                    if (i, index) in Jb_temp:
                        Jb_temp.remove((i, index))
                        no_lines_to_delete = False
        # столбцы
        # ищем в столбце все базисные позиции
        for j in range(n_x):
            basis_indexes_count = 0
            for i in range(m_x):
                if (i, j) in Jb_temp:
                    basis_indexes_count += 1
            # удаляем
            if basis_indexes_count <= 1:
                for index in range(n):
                    if (index, j) in Jb_temp:
                        Jb_temp.remove((index, j))
                        no_rows_to_delete = False
    return Jb_temp


# Теперь расставляем плюсики и минусики.
# Новой базисной позиции (i, j) ставим +.
# А дальше ставим так, чтобы плюсы и минусы чередовались.
# theta выбирается так, чтобы не уйти в отрицательные числа.
def make_plus_minus(i, j, n, m, J_b_temp, x):
    # сюда будем записывать элементы с "+" или с "-"
    plus, minus = [], []
    cycle = 0
    theta = math.inf
    print(x)
    print(J_b_temp)
    # продолжаем до тех пор, пока не останеться элементов в J_b_temp
    while len(J_b_temp) != 0:
        # расставляем + или -
        if cycle % 2 == 0:
            plus.append((i, j))
        else:
            minus.append((i, j))
            if x[i][j] < theta:
                theta = x[i][j]
        J_b_temp.remove((i, j))
        found_vertical = False
        for index in range(n):
            if (index, j) in J_b_temp:
                i = index
                found_vertical = True
                break
        if not found_vertical:
            for index in range(m):
                if (i, index) in J_b_temp:
                    j = index
                    break
        cycle += 1
    return plus, minus, theta


# Обновляем x и J_b
# Там где плюс, мы theta прибавляем, где минус - вычитаем.
def update_x_and_jb(minus, plus, theta, J_b, x):
    is_deleted = False
    for k in minus:
        i, j = k[0], k[1]
        x[i][j] -= theta
        if x[i][j] == 0 and not is_deleted:
            J_b.remove((i, j))
            is_deleted = True
    for k in plus:
        i, j = k[0], k[1]
        x[i][j] += theta
        if k not in J_b:
            J_b.append(k)
    return x, J_b


def matrix_transport_problem(a, b, c):
    # условие баланса
    a, b, c = balance_condition(a, b, c)
    # метод северо-западного угла
    x, J_b = northwest_angle_method(a, b)
    # находим m, n (по новому c, если оно было изменено в условии баланса)
    n, m = np.shape(c)
    # фаза 2
    while True:
        # Нахождение u и v путем решения системы
        u, v = potentials_method(c, J_b)
        # Проверка выполнения условия оптимальности
        i, j = check_optimality_condition(u, v, c, J_b)
        # если условие выполнилось, возвращаем x
        if (i, j) == (None, None):
            return x
        # проверяем угловые ли позиции или не угловые (возвращаем новый Jb с угловыми позициями)
        J_b_temp = find_matrix_corner(m, n, J_b, i, j, x)

        plus, minus, theta = make_plus_minus(i, j, n, m, J_b_temp, x)
        # print('Array with \'pluses\': ' + str(pluses))
        # print('Array with \'minuses\': ' + str(minuses))
        # print('Theta is: ' + str(theta))

        # Добавление или вычитание теты (создание нового X и J_b)
        '''is_deleted = False
        for k in minuses:
            i, j = k[0], k[1]
            x[i][j] -= theta
            if x[i][j] == 0 and not is_deleted:
                J_b.remove((i, j))
                is_deleted = True
        for k in pluses:
            i, j = k[0], k[1]
            x[i][j] += theta
            if k not in J_b:
                J_b.append(k)'''
        x, J_b = update_x_and_jb(minus, plus, theta, J_b, x)
        # print('Промежуточное X:\n' + str(x))
        # print('Промежуточное J basis:\n' + str(J_b))


if __name__ == '__main__':
    # проверка с нулями (важно)
    '''
    a = [0, 0, 0]
    b = [0, 0, 0]
    c = np.array([[0, 0, 0], 
                  [0, 0, 0], 
                  [0, 0, 0]])
    '''
    # var 14
    '''
    a = [100, 70, 170]
    b = [140, 100, 100]
    c = np.array([[1, 5, 6],
                  [2, 4, 2],
                  [4, 1, 1]])
    '''
    # var 11
    a = [50, 170, 80]
    b = [200, 50, 50]
    c = np.array([[2, 5, 6],
                  [3, 4, 3],
                  [6, 1, 2]])
    x = matrix_transport_problem(a, b, c)
    print(x)
