import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def sol_ex(x):
    return np.sin(x) + x + 1

def func(x):
    return np.sin(x)

def get_tridiagonal_matrix(n: int):
    matrix = np.zeros((n, n))
    np.fill_diagonal(matrix, 2)

    for i in range(n - 1):
        matrix[i, i + 1] = -1
        matrix[i + 1, i] = -1
    
    return matrix

def solve_by_finding_inversed_matrix(matrix, vector):
    inversed = np.linalg.inv(matrix)
    return np.dot(inversed, vector)

def solve_by_tridiagonal_matrix(matrix, vector):
    n = len(matrix)
    a = np.zeros(n)
    b = np.zeros(n - 1)
    c = np.zeros(n - 1)

    for i in range(n):
        a[i] = matrix[i, i]
        if i > 0:
            b[i - 1] = matrix[i - 1, i]
            c[i - 1] = matrix[i, i - 1]
    
    return thomas(a, b, c, vector)

def thomas(a, b, c, vector):
    n = len(a)

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

def build_rhs_vector(n, cell_size, a, b):
    vector = np.zeros(n - 1)
    for i in range(n - 1):
        x = cell_size * i
        vector[i] = func(x) * cell_size ** 2
    vector[0] += a
    vector[-1] += b
    return vector

def compute_error_c_norm(solution):
    y_exact = sol_ex(solution)
    return np.max(np.abs(y_exact - solution))

# Параметры задачи
lengh = 10
a = sol_ex(0)
b = sol_ex(lengh)
frames_data = []
errors_inv = []
errors_thomas = []

ns = [2**x for x in range(2, 14)]
for n in ns:
    # решение задачи
    cell_size = lengh / n
    x = np.linspace(0, lengh, n - 1)

    matrix = get_tridiagonal_matrix(n - 1)
    diagonal = np.full(n - 1, 2.)
    upper_diagonal = np.full(n - 2, -1.)
    lower_diagonal = np.full(n - 2, -1.)

    vector = build_rhs_vector(n, cell_size, a, b)
    sol_inv = solve_by_finding_inversed_matrix(matrix, vector)

    vector = build_rhs_vector(n, cell_size, a, b)
    sol_thomas = thomas(diagonal, upper_diagonal, lower_diagonal, vector)

    x_exact = np.linspace(0, lengh, 500)
    y_exact = sol_ex(x_exact)

    frames_data.append((x, sol_inv, sol_thomas, x_exact, y_exact, n))

    # Вычисление ошибки в норме C
    errors_inv.append(compute_error_c_norm(sol_inv))
    errors_thomas.append(compute_error_c_norm(sol_thomas))


# Создание графика решений
fig, ax = plt.subplots(figsize=(10, 6))
line_inv, = ax.plot([], [], label="Inversed Matrix", color="blue", linestyle='--')
line_thomas, = ax.plot([], [], label="Thomas Algorithm", color="green")
line_exact, = ax.plot([], [], label="Exact Solution", color="black", linestyle='-.')

ax.set_xlabel("x")
ax.set_ylabel("Solution")
ax.set_title("Solution of the system")
ax.legend(loc='upper left')
ax.grid(True)

def init():
    ax.set_xlim(-1, 11)
    ax.set_ylim(0, 12)
    return line_inv, line_thomas, line_exact

def update(frame):
    x, sol_inv, sol_thomas, x_exact, y_exact, n_val = frames_data[frame]
    line_inv.set_data(x, sol_inv)
    line_thomas.set_data(x, sol_thomas)
    line_exact.set_data(x_exact, y_exact)
    ax.set_title(f"Solution of the system (n = {n_val})")
    return line_inv, line_thomas, line_exact

ani = FuncAnimation(fig, update, frames=len(frames_data),
                    init_func=init, blit=True, interval=100)

plt.show()

# Лог-график ошибки
plt.figure(figsize=(10, 6))
plt.loglog(ns, errors_inv, label="Inversed Matrix", color="blue", linestyle='--')
plt.loglog(ns, errors_thomas, label="Thomas Algorithm", color="green")
plt.xlabel("Number of cells (n)")
plt.ylabel("Error C-norm")
plt.title("Error C-norm")
plt.legend()
plt.grid(True, which="both", linestyle='--')
plt.show()