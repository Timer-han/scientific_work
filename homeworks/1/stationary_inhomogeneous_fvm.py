import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Пример 1. Всё хорошо
def get_h(x, t):
    return np.cos(x) + x - np.cos(t)
def get_f(x, t):
    return 2 * np.sin(x) * np.cos(x) - np.cos(x) + np.sin(x) * np.sin(t)
def get_k(x):
    return np.sin(x)
def get_s(x):
    return np.sin(x)
def get_h_initial(x):
    return get_h(x, 0)
def get_h_boundaries(t, length):
    return get_h(0, t), get_h(length, t)

# # Пример 2. Решение очень красиво расходится
# def get_h(x, t):
#     return np.sin(x) + x - np.cos(t)
# def get_f(x, t):
#     return np.sin(x) + (x - 0.9) * np.sin(t)
# def get_k(x):
#     return 1
# def get_s(x):
#     return x - 0.9
# def get_h_initial(x):
#     return get_h(x, 0)
# def get_h_boundaries(t, length):
#     return get_h(0, t), get_h(length, t)

# # Пример 3. k = 1, s = 1
# def get_h(x, t):
#     return np.sin(x) + x - np.cos(t)
# def get_f(x, t):
#     return np.sin(x) + np.sin(t)
# def get_k(x):
#     return 1
# def get_s(x):
#     return 1
# def get_h_initial(x):
#     return get_h(x, 0)
# def get_h_boundaries(t, length):
#     return get_h(0, t), get_h(length, t)


def get_tridiagonal_matrix(n: int):
    matrix = np.zeros((n, n))
    np.fill_diagonal(matrix, 2)

    for i in range(n - 1):
        matrix[i, i + 1] = -1
        matrix[i + 1, i] = -1
    
    matrix[0, 0] = 3
    matrix[n - 1, n - 1] = 3
    return matrix

def get_diagonal(n, cell_size, time_interval):
    diagonal = np.full(n - 1, 0.)
    
    for i in range(1, n - 2):
        a = get_k((i - 0.5) * cell_size)
        b = get_k((i + 0.5) * cell_size)
        c = get_k((i + 1.5) * cell_size)
        s = get_s((i + 0.5) * cell_size)
        diagonal[i] = (c * (a + b) + a * (b + c)) / ((a + b) * (b + c)) + s * cell_size ** 2 / (2 * b * time_interval)

    a = get_k(0)
    b = get_k(0.5 * cell_size)
    c = get_k(1.5 * cell_size)
    s = get_s(0.5 * cell_size)
    diagonal[0] = (c * (a + b) + 2 * a * (b + c)) / ((a + b) * (b + c)) + s * cell_size ** 2 / (2 * b * time_interval)

    a = get_k((n - 2.5) * cell_size)
    b = get_k((n - 1.5) * cell_size)
    c = get_k((n - 1.0) * cell_size)
    s = get_s((n - 1.5) * cell_size)
    diagonal[-1] = (2 * c * (a + b) + a * (b + c)) / ((a + b) * (b + c)) + s * cell_size ** 2 / (2 * b * time_interval)

    return diagonal

def get_upper_diagonal(n):
    upper_diagonal = np.full(n - 2, 0.)
    for i in range(0, n - 2):
        a = get_k((i + 0.5) * cell_size)
        b = get_k((i + 1.5) * cell_size)
        upper_diagonal[i] = - b / (a + b)
    return upper_diagonal

def get_lower_diagonal(n):
    lower_diagonal = np.full(n - 2, 0.)
    for i in range(0, n - 2):
        a = get_k((i + 0.5) * cell_size)
        b = get_k((i + 1.5) * cell_size)
        lower_diagonal[i] = - a / (a + b)
    return lower_diagonal

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

def build_rhs_vector(n, cell_size, time_interval, time_step, h_prev):
    vector = np.zeros(n - 1)
    for i in range(0, n - 1):
        x = (i + 0.5) * cell_size
        t = time_step * time_interval
        s = get_s(x)
        k = get_k(x)
        vector[i] = get_f(x, t) * cell_size ** 2 / (2 * k) + s * cell_size ** 2 / (2 * k * time_interval) * h_prev[i]

    lengh = n * cell_size
    h_left, h_right = get_h_boundaries(time_step * time_interval, lengh)
    vector[0] += 2 * h_left * get_k(0) / (get_k(0) + get_k(0.5 * cell_size))
    vector[-1] += 2 * h_right * get_k(lengh) / (get_k(lengh) + get_k((n - 0.5) * cell_size))
    return vector

def compute_error_c_norm(x_exact, solution, time_step, time_interval):
    y_exact = get_h(x_exact, time_step * time_interval)
    return np.max(np.abs(y_exact - solution))

# Параметры задачи
lengh = 10
n = 1000
cell_size = lengh / n
x = np.linspace(cell_size / 2, lengh - cell_size / 2, n - 1)

time_interval = 0.1

frames_data = []
errors_inv = []
errors_thomas = []
h_prev = get_h_initial(x)

for time_step in range(1, 1000):
    # решение задачи

    diagonal = get_diagonal(n, cell_size, time_interval)
    upper_diagonal = get_upper_diagonal(n)
    lower_diagonal = get_lower_diagonal(n)


    # vector = build_rhs_vector(n, cell_size, a, b)
    # sol_inv = solve_by_finding_inversed_matrix(matrix, vector)

    vector = build_rhs_vector(n, cell_size, time_interval, time_step, h_prev)
    sol_thomas = thomas(diagonal, upper_diagonal, lower_diagonal, vector)

    x_exact = np.linspace(cell_size / 2, lengh - cell_size / 2, 1000)
    y_exact = get_h(x_exact, time_step * time_interval)

    # frames_data.append((x, sol_inv, sol_thomas, x_exact, y_exact, n))
    frames_data.append((x, sol_thomas, x_exact, y_exact, n))

    # Вычисление ошибки в норме C
    # errors_inv.append(compute_error_c_norm(x, sol_inv))
    errors_thomas.append(compute_error_c_norm(x, sol_thomas, time_step, time_interval))

    h_prev = sol_thomas


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
    # x, sol_inv, sol_thomas, x_exact, y_exact, n_val = frames_data[frame]
    x, sol_thomas, x_exact, y_exact, n_val = frames_data[frame]
    # line_inv.set_data(x, sol_inv)
    line_thomas.set_data(x, sol_thomas)
    line_exact.set_data(x_exact, y_exact)
    ax.set_title(f"Solution of the system (n = {n_val})")
    return line_inv, line_thomas, line_exact

ani = FuncAnimation(fig, update, frames=len(frames_data),
                    init_func=init, blit=True, interval=100)

plt.show()

# # Лог-график ошибки
# plt.figure(figsize=(10, 6))
# # plt.loglog(ns, errors_thomas, label="Thomas Algorithm", color="red")
# # plt.loglog(ns, errors_inv, label="Inversed Matrix", color="blue", linestyle='--')
# plt.xlabel("Number of cells (n)")
# plt.ylabel("Error C-norm")
# plt.title("Error C-norm")
# plt.legend()
# plt.grid(True, which="both", linestyle='--')
# plt.show()