import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import sys

IC = -8

# Пример 2. Celia
def get_h(x, t):
    return x * 0 + IC
def get_k(x):
    return 7.96 # m/day

# Определить аппроксимацию К в узле
def get_k_node(i, dx, n):
    if i == 0:
        return get_k(dx * 0.5)
    if i == n:
        return get_k(dx * (n-0.5))
    Ki   = get_k(dx * (i+0.5))
    Kim1 = get_k(dx * (i-0.5)) # K_{i-1}
    return 2 * Ki*Kim1 / (Ki + Kim1)
    

def get_sstor(x):
    return 1e-6 # 1/m
def get_h_initial(x):
    return get_h(x, 0)
def get_h_boundaries(t, length):
    return IC, 0.25#get_h(0, t), get_h(length, t)
def get_S(h):
    return get_ө(h) / get_phi() 
def get_Kr(h):
    #val = (h/10)**2
    #print(val)
    #return get_S(h) ** 2
    ө_s = get_phi()
    ө_r = 0.102
    n = 2
    m = 1 - 1 / n
    ө = get_ө(h)
    S_l = (ө - ө_r) / (ө_s - ө_r)
    #print(S_l)
    return S_l ** 2# 0.5 * (1 - (1 - S_l ** (1 / m)) ** m) ** 2
def get_ө(h):
    ө_s = get_phi()
    ө_r = 0.102
    a_vg = 3.35
    n = 2
    m = 1 - 1 / n
    #return ө_s
    if h < 0:
        return ө_r + (ө_s - ө_r) / (1 + abs(a_vg * h) ** n) ** m
    else:
        return ө_s
def get_phi():
    return 0.368


def get_diagonal(n, dx, dt, h_init):
    diagonal = np.full(n - 1, 0.)
    
    # На главной диагонали стоит коэффициент при h_i, он составляется из
    # - пространственной аппроксимации
    # - аппроксимация по времени
    
    # 1. Пространственная аппроксимация
    # Внутренние ячейки
    for i in range(1, n - 2):
        k_l = get_k_node(i,   dx, n)
        k_r = get_k_node(i+1, dx, n)
        Kr_r = get_Kr(max(h_init[i], h_init[i + 1]))
        Kr_l = get_Kr(max(h_init[i - 1], h_init[i]))
        #print(f"Kr: {Kr_l}, {Kr_r}")
        diagonal[i] = (k_l * Kr_l + k_r * Kr_r)
        
    # Левая граничная ячейка с номером 0
    k_l = get_k(0.5 * dx) # Берем из ячейки
    k_r = get_k_node(1, dx, n)
    Kr_l = get_Kr(max(h_init[0],-10)) # захардкодил ГУ слева
    Kr_r = get_Kr(max(h_init[0], h_init[1]))
    diagonal[0] = 0.5 * k_l * Kr_l + k_r * Kr_r 

    # Правая граничная ячейка с номером n-1
    k_l = get_k_node(n-1, dx, n)
    k_r = get_k((n-1.5) * dx) # из ячейки
    Kr_l = get_Kr(max(h_init[-2], h_init[-1]))
    Kr_r = get_Kr(max(h_init[-1], 0.25)) # захардкодил ГУ справа
    diagonal[-1] = k_l * Kr_l + 0.5 * k_r * Kr_r 

    # 2. Aппроксимация по времени
    for i in range(0,n-1):
        S            = get_S(h_init[i])
        sstor        = get_sstor((i + 0.5) * dx)
        diagonal[i] += (dx ** 2) * sstor * S / dt;
    
    #print(diagonal)
    #exit(-1)

    return diagonal

def build_rhs_vector(n, dx, dt, time_step, h_prev, h_init):
    vector = np.zeros(n)

    for i in range(0, n - 1):
        x = (i + 0.5) * dx
        t = time_step * dt
        sstor = get_sstor(x)
        k = get_k(x)
        S = get_S(h_init[i])
        vector[i] = sstor * S * (dx**2) / dt * h_prev[i]

    length = n * dx
    h_left, h_right = get_h_boundaries(time_step * dt, length)
    
    # Левая граничная ячейка с номером 0
    k_l = get_k(0.5 * dx) # Берем из ячейки
    Kr_l = get_Kr(max(h_init[0],-10)) # захардкодил ГУ слева
    vector[0] += 0.5 * k_l * Kr_l * h_left 

    # Правая граничная ячейка с номером n-1
    k_r = get_k((n-1.5) * dx) # из ячейки
    Kr_r = get_Kr(max(h_init[-1], 0.25)) # захардкодил ГУ справа
    vector[-1] += 0.5 * k_r * Kr_r * h_right 
    
    #vector[0] += 2 * h_left * get_k(0) / (get_k(0) + get_k(0.5 * dx))
    #vector[-1] += 2 * h_right * get_k(length) / (get_k(length) + get_k((n - 0.5) * dx))
    
    return vector

def get_upper_diagonal(n, h_init):
    upper_diagonal = np.full(n - 2, 0.)
    for i in range(0, n - 2):
        a = get_k((i + 0.5) * dx)
        b = get_k((i + 1.5) * dx)
        k = get_k_node(i+1, dx, n)
        Kr = get_Kr(max(h_init[i], h_init[i + 1]))
        upper_diagonal[i] = -k * Kr#- b / (a + b) * Kr
    #print(upper_diagonal)
    #exit(-1)
    return upper_diagonal

def get_lower_diagonal(n, h_init):
    lower_diagonal = np.full(n - 2, 0.)
    for i in range(0, n - 2):
        a = get_k((i + 0.5) * dx)
        b = get_k((i + 1.5) * dx)
        k = get_k_node(i, dx, n)
        Kr = get_Kr(max(h_init[i], h_init[i + 1]))
        lower_diagonal[i] = -k * Kr#- a / (a + b) * Kr
    #print(lower_diagonal)
    #exit(-1)
    return lower_diagonal

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

# Функция вычисления нормы невязки r = ||b(h) - A(h)*h||_2
# Аргументы: центральная диагональ, нижняя диагональ, вектор правой части (right-hand side), решение, dt
def compute_residual_norm(n, diag_c, diag_u, diag_l, rhs, h, dt):
    #print(diag_c)
    #print(diag_u)
    #print(diag_l)
    #print(h)
    #exit(1)
    r_norm = 0
    for i in range(0, n-2):
        r_i = rhs[i] - diag_c[i] * h[i] # главная диагональ
        if i > 0: 
            r_i -= diag_l[i-1] * h[i-1] # нижняя диагональ
        if i < n-2: 
            r_i -= diag_u[i]   * h[i+1] # верхняя диагональ
        #print(r_i)
        r_norm += r_i**2
        #print(f"r[{i}] = {r_i}")
    return r_norm ** 0.5
            
    

# Параметры задачи
eps_abs = 1e-5
eps_rel = 1e-5
length = 1 # m
n = 32
dx = length / n
x = np.linspace(dx / 2, length - dx / 2, n - 1)
#print(x)

dt = 1e-9 # days
T  = 1e-7
nt = int(T/dt)
save_intensity = int(nt/10)

frames_data = []
errors_inv = []
errors_thomas = []

# Начальные условия - напор в начальный момент
h_init = get_h_initial(x)

# Решение на предыдущем шаге по времени - h^n
# Для понятности переименовал его из h_prev в h_n
h_n = h_init.copy()

# Решение на новом шаге по времени - h^{n+1}
# Для понятности назвал его просто h
h = h_n.copy()

for time_step in range(1, nt):
    # решение задачи
    print(f"-- Time step {time_step}, t_beg = {(time_step-1) * dt}")

    converged = False
    r_0 = 0
    # Цикл метода простой итерации (Picard method)
    for k in range(100):
        # Формируем матрицу, записывая диагонали
        diagonal = get_diagonal(n, dx, dt, h)
        upper_diagonal = get_upper_diagonal(n, h)
        lower_diagonal = get_lower_diagonal(n, h)
        #print(diagonal)
        #print(lower_diagonal)
        #print(upper_diagonal)

        # Составляем вектор правой части
        vector = build_rhs_vector(n, dx, dt, time_step, h_n, h)
        
        # Считаем невязку
        # В данном случае h выступает как h^{n+1,k}
        # Переход к h^{n+1,k+1} случится после вызова прогонки
        r_k = compute_residual_norm(n, diagonal, upper_diagonal, lower_diagonal, vector, h, dt)
        if k == 0:
            r_0 = r_k
        print(f"  iter {k}, r = {r_k}")
        
        # Смотрим на сходимость 
        if (r_k < eps_abs or r_k < eps_rel * r_0) and k > 0:
            print("  Converged")
            converged = True
            break
        
        # Решение системы прогонкой
        # h становится h^{n+1,k+1}
        h = thomas(diagonal, upper_diagonal, lower_diagonal, vector)
        #print(h)
        #exit(1)

    if not converged:
        print(f"Не сошлось на шаге {time_step}!")
        break
    
    # Сошлись, обновляем решение на предыдущем шаге
    h_n = h
    #print(h)

    x_exact = np.linspace(dx / 2, length - dx / 2, 1000)
    y_exact = get_h(x_exact, time_step * dt)

    # frames_data.append((x, sol_inv, sol_thomas, x_exact, y_exact, n))
    sat = h.copy()
    for i in range(0,n-1):
        sat[i] = get_S(h[i])
    
    if time_step % save_intensity == 1:
        frames_data.append((x, h, h_init, x_exact, y_exact, n))
    #frames_data.append((x, h, n))

    # Вычисление ошибки в норме C
    # errors_inv.append(compute_error_c_norm(x, sol_inv))
    # errors_thomas.append(compute_error_c_norm(x, sol_thomas, time_step, dt))

    h_init = h
print(f"save_intensity = {save_intensity}")


# Создание графика решений
fig, ax = plt.subplots(figsize=(10, 6))
line_inv, = ax.plot([], [], label="Previous step solution", color="blue", linestyle='--')
line_thomas, = ax.plot([], [], label="Thomas Algorithm", color="green")
line_exact, = ax.plot([], [], label="Initial condition", color="black", linestyle='-.')

ax.set_xlabel("x")
ax.set_ylabel("Solution")
ax.set_title("Solution of the system")
ax.legend(loc='upper left')
ax.grid(True)




def init():
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(1.1*IC, 1)
    return line_inv, line_thomas, line_exact

def update(frame):
    # x, sol_inv, sol_thomas, x_exact, y_exact, n_val = frames_data[frame]
    x, sol_thomas, h_prev, x_exact, y_exact, n_val = frames_data[frame]
    # line_inv.set_data(x, sol_inv)
    line_inv.set_data(x, h_prev)
    line_thomas.set_data(x, sol_thomas)
    line_exact.set_data(x_exact, y_exact)
    ax.set_title(f"Solution of the system (n = {n_val})")
    return line_inv, line_thomas, line_exact




ani = FuncAnimation(fig, update, frames=len(frames_data),
                    init_func=init, blit=True, interval=100)

# pillow_writer = PillowWriter(fps=10)
# ani.save('animation.gif', writer=pillow_writer)
plt.show()
