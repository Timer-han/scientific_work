import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import sys
import os

IC = -8

# Пример 2. Celia
def get_h(x, t):
    return x * 0 + IC
def get_k(x):
    return 7.96 # m/day

# Пока что зависимости от координаты нет,
# то есть область однородная
def get_VG_param(x):
    th_s = get_phi()
    th_r = 0.102
    n = 2
    m = 1 - 1 / n
    a_vg = 3.35*10
    return th_s, th_r, n, m, a_vg
    

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
    ө_s, ө_r, n, m, a = get_VG_param(0)
    ө = get_ө(h)
    S_l = (ө - ө_r) / (ө_s - ө_r)
    #print(S_l)
    return S_l ** 2# 0.5 * (1 - (1 - S_l ** (1 / m)) ** m) ** 2
def get_ө(h):
    ө_s, ө_r, n, m, a_vg = get_VG_param(0)
    if h < 0:
        return ө_r + (ө_s - ө_r) / (1 + abs(a_vg * h) ** n) ** m
    else:
        return ө_s
    
def get_theta_der(h):
    theta_s, theta_r, n, m, alpha = get_VG_param(0)
    if h < 0.0:
         return (theta_s - theta_r) * m*n * alpha * pow(-alpha*h, n-1.0) * pow(1.0 + pow(-alpha*h, n), -m - 1.0);
    else:
         return 0.0;
    
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
        coef         = S * sstor
        coef        += get_theta_der(h_init[i])
        diagonal[i] += coef * (dx ** 2) / dt
    
    #print(diagonal)
    #exit(-1)

    return diagonal

def build_rhs_vector(n, dx, dt, time_step, h_prev, h_init):
    vector = np.zeros(n)

    # Часть, связанная с производной по времени
    for i in range(0, n - 1):
        x = (i + 0.5) * dx
        t = time_step * dt
        sstor = get_sstor(x)
        k = get_k(x)
        S = get_S(h_init[i])
        vector[i] = sstor * S * (dx**2) / dt * h_prev[i]
        theta_n = get_ө(h_prev[i])
        theta_k = get_ө(h_init[i])
        c_theta = get_theta_der(h_init[i])
        vector[i] += c_theta * h_init[i] * (dx**2) / dt 
        vector[i] -= (theta_k - theta_n) * (dx**2) / dt 
        
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

# --------------------------------------------------------------------
# Инициализация некоторых переменных
frame_top = 1
frame_bottom = 1.1*IC

# --------------------------------------------------------------------
# Чтение данных из файлов или их вычисление

def read_data(filename):
    data = np.load(filename, allow_pickle=True)
    # Get the number of tuples from the file names
    max_i = max(int(key.split('_')[1]) for key in data.files)
    # Get the number of elements in each tuple
    max_j = max(int(key.split('_')[2]) for key in data.files)

    print(f"Reading data from {filename}, max_i = {max_i}, max_j = {max_j}")
    
    result = []
    for i in range(max_i + 1):
        tuple_data = []
        for j in range(max_j + 1):
            key = f'data_{i}_{j}'
            if key in data:
                tuple_data.append(data[key])
        result.append(tuple(tuple_data))
        print(f"Tuple {i}: {tuple_data[:3]}")
    print(f"Readed {len(result)} tuples from {filename}")
    return result

answer = False
if (os.path.exists("frames_data.npz")
    and os.path.exists("low_eps_solution_compare.npz")
    and os.path.exists("high_eps_solution_compare.npz")):
    print("Found existing data files, load them?")
    while True:
        answer = input("Type 'yes' to load or 'no' to compute new data: ").strip().lower()
        if answer in ['yes', 'no', 'y', 'n']:
            if answer in ['yes', 'y']:
                print("Loading existing data...")
                answer = True
                break
            else:
                print("Computing new data...")
                answer = False
                break
        print("Please type 'yes' or 'no'.")
    

if answer:
    frames_data = read_data("frames_data.npz")
    low_eps_solution_compare = read_data("low_eps_solution_compare.npz")
    high_eps_solution_compare = read_data("high_eps_solution_compare.npz")
    
else:

    # Параметры задачи

    length = 1 # m
    T  = 1e-2

    n_array = [2**i for i in range(4, 13)] # Количество ячеек
    t_array = [1e-2 / 10 / 2**i for i in range(1, 7)] # Временные шаги
    t_array.extend([t_array[-1], t_array[-1], t_array[-1]])

    if len(n_array) != len(t_array):
        print("Количество ячеек и временных шагов не совпадает!")
        print(f"n_array: {len(n_array)}, t_array: {len(t_array)}")
        sys.exit(1)

    

    frames_data = []
    errors_inv = []
    errors_thomas = []
    sum_flows = []

    low_eps_solution_compare = []
    high_eps_solution_compare = []

    for i in range(len(n_array)):
        n = n_array[i]
        dx = length / n
        dt = t_array[i]
        nt = int(T/dt)
        save_intensity = int(nt/20)
        x = np.linspace(dx / 2, length - dx / 2, n - 1)

        for eps in [1e-1, 1e-7]:
            eps_abs = eps
            eps_rel = eps

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
                
                x_exact = np.linspace(dx / 2, length - dx / 2, 1000)
                y_exact = get_h(x_exact, time_step * dt)

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

                # frames_data.append((x, sol_inv, sol_thomas, x_exact, y_exact, n))
                sat = h.copy()
                for i in range(0,n-1):
                    sat[i] = get_S(h[i])


                flow_in_cells = []
                a = - get_k(dx * 0) * get_k_node (0, dx, n) * (h[1] - h[0]) / dx
                for i in range(1, n):
                    b = - get_k(dx * i) * get_k_node (0, dx, n) * (h[1] - h[0]) / dx
                    flow_in_cells.append(a + b)
                    a = b
                flow_in_cells = flow_in_cells - np.average(flow_in_cells)


                if time_step % save_intensity == 1:
                    frames_data.append((x, h, h_init, x_exact, y_exact, flow_in_cells, n))
                    frame_bottom = min(np.min(h), frame_bottom)
                    frame_top = max(np.max(flow_in_cells + h), frame_top)
                #frames_data.append((x, h, n))

                # Вычисление ошибки в норме C
                # errors_inv.append(compute_error_c_norm(x, sol_inv))
                # errors_thomas.append(compute_error_c_norm(x, sol_thomas, time_step, dt))



                h_init = h

            if eps == 1e-1:
                high_eps_solution_compare.append((x, h, x_exact, y_exact, n))
            else:
                low_eps_solution_compare.append((x, h, x_exact, y_exact, n))
# print(f"save_intensity = {save_intensity}")


for i in frames_data:
    print(f"flow difference on step i: {np.average(i[5])}")



# --------------------------------------------------------------------
# Создание графика решений
fig, ax = plt.subplots(figsize=(10, 6))
line_prev, = ax.plot([], [], label="Previous step solution", color="blue", linestyle='--')
line_thomas, = ax.plot([], [], label="Thomas Algorithm", color="green")
line_exact, = ax.plot([], [], label="Exact solution", color="black", linestyle='-.')
line_flow, = ax.plot([], [], label="Flow difference in cells", color="red", linestyle="-")

ax.set_xlabel("x")
ax.set_ylabel("Solution")
ax.set_title("Solution of the system")
ax.legend(loc='upper left')
ax.grid(True)



def init():
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(frame_bottom * 1.1, frame_top * 1.1)
    return line_prev, line_thomas, line_exact

def update(frame):
    # x, sol_inv, sol_thomas, x_exact, y_exact, n_val = frames_data[frame]
    x, sol_thomas, h_prev, x_exact, y_exact, flow_in_cells, n_val = frames_data[frame]
    # line_prev.set_data(x, sol_inv)
    line_prev.set_data(x, h_prev)
    line_thomas.set_data(x, sol_thomas)
    line_exact.set_data(x_exact, y_exact)
    line_flow.set_data(x, flow_in_cells)
    ax.set_title(f"Solution of the system (n = {n_val})")
    return line_prev, line_thomas, line_exact

ani = FuncAnimation(fig, update, frames=len(frames_data),
                    init_func=init, blit=True, interval=100)

pillow_writer = PillowWriter(fps=10)

print("Saving animation.gif...")
ani.save('animation.gif', writer=pillow_writer)
plt.show()


# --------------------------------------------------------------------
# Сравнение решений с разными eps
fig, ax = plt.subplots(figsize=(10, 6))
line_exact, = ax.plot([], [], label="Initial condition", color="black", linestyle='-.')
line_low_eps_solution, = ax.plot([], [], label="Low eps solution", color="orange", linestyle=':')
line_high_eps_solution, = ax.plot([], [], label="High eps solution", color="purple", linestyle=':')

ax.set_xlabel("x")
ax.set_ylabel("Solution")
ax.set_title("Solution of the system")
ax.legend(loc='upper left')
ax.grid(True)


def init_sol_comp():
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(frame_bottom * 1.1, frame_top * 1.1)

    return line_high_eps_solution, line_low_eps_solution, line_exact, ax.title

def update_sol_comp(frame):
    x, h, x_exact, y_exact, n_val = high_eps_solution_compare[frame]
    line_high_eps_solution.set_data(x, h)

    x, h, x_exact, y_exact, n_val = low_eps_solution_compare[frame]
    line_low_eps_solution.set_data(x, h)

    line_exact.set_data(x_exact, y_exact)

    ax.set_title(f"Solution of the system (n = {n_val})")
    
    return line_high_eps_solution, line_low_eps_solution, line_exact, ax.title


different_eps_compare_ani = FuncAnimation(fig, update_sol_comp, frames=len(high_eps_solution_compare),
                    init_func=init_sol_comp, blit=True, interval=500)

pillow_writer = PillowWriter(fps=2)
print("Saving different_eps_compare.gif...")
different_eps_compare_ani.save('different_eps_compare.gif', writer=pillow_writer)

plt.show()


# --------------------------------------------------------------------
# Вычисление разницы между решениями с разными eps
solution_diff = []
for i in range(len(high_eps_solution_compare)):
    x, h_high, x_exact, y_exact, n_val = high_eps_solution_compare[i]
    x, h_low, _, _, _ = low_eps_solution_compare[i]
    solution_diff.append((n_val, np.max(np.abs(h_high - h_low))))

# Log-log plot of solution differences
plt.figure(figsize=(12, 6))
plt.loglog([n for n, diff in solution_diff], [diff for n, diff in solution_diff])
plt.xlabel("Number of cells (n)")
plt.ylabel("Max difference between solutions")
plt.title("Max difference between high and low eps solutions")
plt.grid(True)
print("Saving solution_difference_loglog.png...")
plt.savefig("solution_difference_loglog.png")
plt.show()


# --------------------------------------------------------------------
# Сравнение решений с разными n и eps на одном графике
plt.figure(figsize=(10, 6))

for i in range(len(high_eps_solution_compare)):
    x, h_low, _, _, _ = low_eps_solution_compare[i]
    plt.plot(x, h_low, label=f"Low eps solution (n={n_val})", linestyle='-', color=plt.cm.hsv(i / len(low_eps_solution_compare)))

    x, h_high, x_exact, y_exact, n_val = high_eps_solution_compare[i]
    plt.plot(x, h_high, label=f"High eps solution (n={n_val})", linestyle=':', color=plt.cm.hsv(i / len(low_eps_solution_compare) + 0.03))

plt.xlabel("x")
plt.ylabel("Solution")
plt.title("Comparison of solutions with different n and eps")
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.grid(True)
plt.tight_layout()
print("Saving comparison_n_eps_solutions.png...")
plt.savefig("comparison_n_eps_solutions.png", bbox_inches='tight')
plt.show()

# --------------------------------------------------------------------
# Сохранение данных

def save_data(filename, data_list):
    data_dict = {}
    for i, tup in enumerate(data_list):
        for j, arr in enumerate(tup):
            data_dict[f'data_{i}_{j}'] = np.array(arr)
    np.savez(filename, **data_dict)

# Save all data
if not answer:
    print("Saving data to files...")
    save_data("frames_data.npz", frames_data)
    save_data("low_eps_solution_compare.npz", low_eps_solution_compare)
    save_data("high_eps_solution_compare.npz", high_eps_solution_compare)

print("All data saved successfully.")
print("Done!")
