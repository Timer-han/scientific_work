import numpy as np
import matplotlib.pyplot as plt


def func(x):
    return (x - 5)**2

def get_tridiagonal_matrix(n: int):
    matrix = np.zeros((n, n))
    np.fill_diagonal(matrix, 2)

    for i in range(n - 1):
        matrix[i, i + 1] = -1
        matrix[i + 1, i] = -1

    matrix[0, 0] = 3
    matrix[-1, -1] = 3

    return matrix

def solve_by_finding_inversed_matrix(matrix, vector):
    inversed = np.linalg.inv(matrix)
    return np.dot(inversed, vector)

def solve_by_tridiagonal_matrix_algorithm(matrix, vector):
    n = len(matrix)
    a = np.zeros(n)
    b = np.zeros(n - 1)
    c = np.zeros(n - 1)

    for i in range(n):
        a[i] = matrix[i, i]
        if i > 0:
            b[i - 1] = matrix[i - 1, i]
            c[i - 1] = matrix[i, i - 1]
    
    b[0] /= a[0]
    for i in range(1, n - 1):
        b[i] /= (a[i] - c[i - 1] * b[i - 1])
    
    vector[0] /= a[0]
    for i in range(1, n):
        vector[i] = (vector[i] - c[i - 1] * vector[i - 1]) / (a[i] - c[i - 1] * b[i - 1])
    
    x = np.zeros(n)
    x[-1] = vector[-1]
    for i in range(n - 2, -1, -1):
        x[i] = vector[i] - b[i] * x[i + 1]
    
    return x

def residual(matrix, vector, solution):
    return np.dot(matrix, solution) - vector

def get_norm(vector):
    return np.linalg.norm(vector)



def get_initial_data(n, cell_size, a, b):
    vector = np.zeros(n - 1)
    for i in range(n - 1):
        x = cell_size * (i + 0.5)
        vector[i] = func(x) * cell_size ** 2
    vector[0] += 2 * a
    vector[-1] += 2 * b
    return vector

lengh = 10
n = 1000
a = 0
b = 0
cell_size = lengh / n
matrix = get_tridiagonal_matrix(n - 1)

vector = get_initial_data(n, cell_size, a, b)

solution_by_inversed_matrix = solve_by_finding_inversed_matrix(matrix, vector)
solution_by_tridiagonal_matrix_algorithm = solve_by_tridiagonal_matrix_algorithm(matrix, vector)
residual_by_inversed_matrix = residual(matrix, vector, solution_by_inversed_matrix)
residual_by_tridiagonal_matrix_algorithm = residual(matrix, vector, solution_by_tridiagonal_matrix_algorithm)
norm_by_inversed_matrix = get_norm(residual_by_inversed_matrix)
norm_by_tridiagonal_matrix_algorithm = get_norm(residual_by_tridiagonal_matrix_algorithm)
print("Solution by inversed matrix:")
print(solution_by_inversed_matrix[:n:n//10])
print("Norm of residual:")
print(norm_by_inversed_matrix)
print("-" * 50)
print("Solution by tridiagonal matrix algorithm:")
print(solution_by_tridiagonal_matrix_algorithm[:n:n//10])
print("Norm of residual:")
print(norm_by_tridiagonal_matrix_algorithm)
print("-" * 50)

plt.figure(figsize=(10, 6))
plt.plot(np.linspace(0, lengh, n - 1), solution_by_inversed_matrix, label="Inversed Matrix", color="blue")
plt.plot(np.linspace(0, lengh, n - 1), solution_by_tridiagonal_matrix_algorithm, label="Tridiagonal Matrix Algorithm", color="red")
plt.xlabel("x")  # Подпись оси X
plt.ylabel("Solution")  # Подпись оси Y
plt.title("Solution of the system")  # Заголовок графика
plt.legend()  # Показываем легенду
plt.grid(True)  # Включаем сетку
plt.show()  # Показываем график