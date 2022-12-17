import numpy as np
import copy
from collections.abc import Sequence

epsilon = 0.5e-4

matrix_full = [
    [-0.14, 1.21, 1.58, 4.93],
    [0.94, 0.12, -0.61, 1.63],
    [-0.26, -0.55, -0.24, -2.41]
]

matrix_A = [
    [0.94, 0.12, -0.61],
    [-0.26, -0.55, -0.24],
    [-0.14, 1.21, 1.58]
]

vector_b = [
    1.63,
    -2.41,
    4.93
]


# Метод Гаусса
def method_gauss(matrix):
    len_matrix = len(matrix)
    matrix_copy = copy.deepcopy(matrix)

    for i in range(len_matrix):
        for j in range(len_matrix + 1):
            matrix_copy[i][j] = matrix_copy[i][j] / matrix[i][i]
        for j in range(i + 1, len_matrix):
            str_j_div_str_i = matrix_copy[j][i] / matrix_copy[i][i]  # множитель для зануления
            for k in range(len_matrix + 1):
                matrix_copy[j][k] -= matrix_copy[i][k] * str_j_div_str_i
        for j in range(len_matrix):
            for k in range(len_matrix + 1):
                matrix[j][k] = matrix_copy[j][k]

    for i in range(len_matrix - 1, -1, -1):
        for j in range(len_matrix, -1, -1):
            matrix_copy[i][j] = matrix_copy[i][j] / matrix[i][i]
        for j in range(i - 1, -1, -1):
            str_j_div_str_i = matrix_copy[j][i] / matrix_copy[i][i]
            for k in range(len_matrix, -1, -1):
                matrix_copy[j][k] -= matrix_copy[i][k] * str_j_div_str_i

    answer = []
    for i in range(len_matrix):
        answer.append(matrix_copy[i][len_matrix])
    return 'Метод Гаусса:' + '\n' + \
           'x1 = ' + str(answer[0]) + '\n' + \
           'x2 = ' + str(answer[1]) + '\n' + \
           'x3 = ' + str(answer[2]) + '\n'


print(method_gauss(matrix_full))


# Метод Гаусса с выбором главного элемента по всей матрице
def change_order_of_lines(matrix, col):
    max_element = matrix[col][col]
    max_row = col
    for i in range(col + 1, len(matrix)):
        if abs(matrix[i][col]) > abs(max_element):
            max_element = matrix[i][col]
            max_row = i
    if max_row != col:
        matrix[col], matrix[max_row] = matrix[max_row], matrix[col]


def is_only_one_answer(matrix):
    for i in range(len(matrix)):
        if not matrix[i][i]:
            return True
    return False


def method_gauss_with_choice(matrix):
    n = len(matrix)
    for k in range(n - 1):
        change_order_of_lines(matrix, k)
        for i in range(k + 1, n):
            str_i_div_str_k = matrix[i][k] / matrix[k][k]
            matrix[i][-1] -= matrix[k][-1] * str_i_div_str_k
            for j in range(k, n):
                matrix[i][j] -= matrix[k][j] * str_i_div_str_k

    if is_only_one_answer(matrix):
        print('Система имеет бесконечное множество решений.')
        return

    answer = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        answer[i] = (matrix[i][-1] - sum([matrix[i][j] * answer[j] for j in range(i + 1, n)])) / matrix[i][i]

    return 'Метод Гаусса с выбором:' + '\n' + \
           'x1 = ' + str(answer[0]) + '\n' + \
           'x2 = ' + str(answer[1]) + '\n' + \
           'x3 = ' + str(answer[2]) + '\n'


print(method_gauss_with_choice(matrix_full))


# Метод Якоби
def summa(a: Sequence[float], x: Sequence[float], j: int) -> float:
    summa: float = 0
    for i, (m, y) in enumerate(zip(a, x)):
        if i != j:
            summa += m * y
    return summa


def method_jacobi(a: Sequence[Sequence[float]],
                  b: Sequence[float]):
    y = [c / a[i][i] for i, c in enumerate(b)]

    answer = [-(summa(c, y, i) - b[i]) / a[i][i] for i, c in enumerate(a)]
    count_iter = 0

    while max((abs(x0 - y0) for x0, y0 in zip(answer, y))) > epsilon:
        count_iter += 1
        for i, j in enumerate(answer): y[i] = j
        for i, c in enumerate(a): answer[i] = -(summa(c, y, i) - b[i]) / a[i][i]
    return 'Метод Якоби:' + '\n' + \
           'x1 = ' + str(answer[0]) + '\n' + \
           'x2 = ' + str(answer[1]) + '\n' + \
           'x3 = ' + str(answer[2]) + '\n' + \
           'количество итераций: ' + str(count_iter) + '\n'


print(method_jacobi(matrix_A, vector_b))


# Метод Гаусса-Зейделя
def method_gauss_zeidel(A, b):
    n = len(A)
    answer = np.zeros(n)  # нулевой вектор
    count_iter = 0

    converge = False
    while not converge:
        x_n = np.copy(answer)
        for count_iter in range(n):
            s1 = sum(A[count_iter][j] * x_n[j] for j in range(count_iter))
            s2 = sum(A[count_iter][j] * answer[j] for j in range(count_iter + 1, n))
            x_n[count_iter] = (b[count_iter] - s1 - s2) / A[count_iter][count_iter]

        converge = np.linalg.norm(x_n - answer) <= epsilon
        answer = x_n
        count_iter += 1

    return 'Метод Гаусса-Зейделя:' + '\n' + \
           'x1 = ' + str(answer[0]) + '\n' + \
           'x2 = ' + str(answer[1]) + '\n' + \
           'x3 = ' + str(answer[2]) + '\n' + \
           'количество итераций: ' + str(count_iter) + '\n'


print(method_gauss_zeidel(matrix_A, vector_b))
