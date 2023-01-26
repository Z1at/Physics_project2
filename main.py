import numpy as np
import matplotlib.pyplot as plt

dt = 0.001  # step of integration
BETA1 = 0.0001
BETA2 = 3
ALPHA = 10
v = 0.5  # frequency
J0 = 0  # starting magnetization
m0 = 4 * np.pi * 10 ** (-7)  # magnetic constant


# Функция для подсчёта намагниченности на каждом шаге :
def calculate_Jx(Jx, Hx):
    Jxi = Jx[-1] + dt * ALPHA * (Hx[-1] - BETA1 * Jx[-1] - BETA2 * (Jx[-1]) ** 3)
    return Jxi


# Функция для подсчёта магнитной индукции на каждом шаге:
def calculate_Bx(Jx, Hx):
    Bxi = m0 * (Jx[-1] + Hx[-1])
    return Bxi


# Функция для подсчёта напряженности магнитного поля на каждом шаге:
def calculate_Hx(v, t):
    Hxi = np.sin(v * t)
    return Hxi


# Функция для построения петли гистерезиса и кривой первоначальной намагниченности:
def dependence_BH(v):
    t = 0
    Jx = np.array([J0])
    Hx = np.array([np.sin(v * t)])
    Bx = np.array([0])
    while (Hx[-1] <= 1):
        t += dt
        Hxi = calculate_Hx(v, t)
        if (Hx[-1] < Hxi):
            Hx = np.append(Hx, Hxi)
            Jxi = calculate_Jx(Jx, Hx)
            Jx = np.append(Jx, Jxi)
            Bxi = calculate_Bx(Jx, Hx)
            Bx = np.append(Bx, Bxi)
        else:
            break
    while (Hx[-1] >= -1):
        t += dt
        Hxi = calculate_Hx(v, t)
        if (Hx[-1] > Hxi):
            Hx = np.append(Hx, Hxi)
            Jxi = calculate_Jx(Jx, Hx)
            Jx = np.append(Jx, Jxi)
            Bxi = calculate_Bx(Jx, Hx)
            Bx = np.append(Bx, Bxi)
        else:
            break
    while (Hx[-1] <= 1):
        t += dt
        Hxi = calculate_Hx(v, t)
        if (Hx[-1] < Hxi):
            Hx = np.append(Hx, Hxi)
            Jxi = calculate_Jx(Jx, Hx)
            Jx = np.append(Jx, Jxi)
            Bxi = calculate_Bx(Jx, Hx)
            Bx = np.append(Bx, Bxi)
        else:
            break
    plt.plot(Hx, Bx)


# Функция для построения графика зависимости u(H) магнитной проницаемости от напряженности магнитного поля:
def dependence_uH():
    t = 0
    Jx = np.array([J0])
    Hx = np.array([np.sin(v * t)])
    u = np.array([1])
    while (Hx[-1] <= 1):
        t += dt
        Hxi = calculate_Hx(v, t)
        if (Hx[-1] < Hxi):
            Hx = np.append(Hx, Hxi)
            Jxi = calculate_Jx(Jx, Hx)
            Jx = np.append(Jx, Jxi)
            ui = 1 + (Jx[-1] / Hx[-1])
            u = np.append(u, ui)
        else:
            break
    while (Hx[-1] >= -1):
        t += dt
        Hxi = calculate_Hx(v, t)
        if (Hx[-1] > Hxi):
            Hx = np.append(Hx, Hxi)
            Jxi = calculate_Jx(Jx, Hx)
            Jx = np.append(Jx, Jxi)
            ui = 1 + (Jx[-1] / Hx[-1])
            u = np.append(u, ui)
        else:
            break
    k = np.argmax(u)
    H_res = np.array([0])
    u_res = np.array([u[0]])
    for i in range(k, k - 20, -1):
        H_res = np.append(H_res, Hx[i])
        u_res = np.append(u_res, u[i])
    plt.plot(H_res, u_res)


# Функция для построения графика зависимости J(H) намагниченности от напряженности магнитного поля:
def dependence_JH(F, J0, label):
    t = 0
    Jx = np.array([J0])
    Hx = np.array([np.sin(v * t)])
    while (Hx[-1] <= 1):
        t += dt
        Hxi = calculate_Hx(v, t)
        if (Hx[-1] < Hxi):
            Hx = np.append(Hx, Hxi)
            Jxi = calculate_Jx(Jx, Hx)
            Jx = np.append(Jx, Jxi)
        else:
            break
    plt.plot(Hx, Jx, label=label)


# Кривая первоначальной намагниченности и петля гистерезиса
dependence_BH(v)
plt.title(f"Кривая первоначальной намагниченности")
plt.xlabel("H, A/м")
plt.ylabel("B, Тл")
plt.grid()
plt.show()


# Зависимость магнитной проницаемости от напряженности магнитного поля
dependence_uH()
plt.title( f"Зависимость μ(H) магнитной проницаемости \n от напряженности магнитного поля")
plt.xlabel("H, A/м")
plt.ylabel("μ, Гн/м")
plt.grid()
plt.show()


# Зависимость намагниченности от напряженности магнитного поля для различных значений частоты
v = [6, 9, 12]
for v in v:
    dependence_JH(v, J0, f"v={round(v, 2)} Гц")
plt.title( f"Завивисимость J(H) намагниченности \n от напряженности внешнего поля при различных частотах")
plt.xlabel("H, A/м")
plt.ylabel("J, А/м")
plt.legend()
plt.grid()
plt.show()


# Зависимость намагниченности от напряженности магнитного поля для различных значений стартовой намагниченности
J0 = [0, 0.05, 0.1, 0.3]
v = 12
for J0 in J0:
    dependence_JH(v, J0, f"J0={J0} А/м")
plt.title( f"Завивисимость J(H) намагниченности от напряженности \n внешнего поля при различных начальных намагничиваниях")
plt.xlabel("H, A/м")
plt.ylabel("J, А/м")
plt.legend()
plt.grid()
plt.show()
