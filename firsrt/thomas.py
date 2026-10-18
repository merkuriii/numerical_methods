"""Лабораторная 1.2. Метод прогонки (обозначения методички Пирумова, §1.5, формулы (1.6)-(1.8)).
Матрица целиком НЕ хранится: только три диагонали.

Во входном файле система записана как в задании:  a_i*x_{i-1} + b_i*x_i + c_i*x_{i+1} = d_i,  a_1 = c_n = 0.
Внутри программа переходит к канонической форме методички:  a_i*x_{i-1} - B_i*x_i + c_i*x_{i+1} = d_i,  B_i = -b_i,
x_i = P_{i+1}*x_{i+1} + Q_{i+1},  P_1 = Q_1 = 0,
P_{i+1} = c_i / (B_i - a_i*P_i),   Q_{i+1} = (a_i*Q_i - d_i) / (B_i - a_i*P_i),   x_n = Q_{n+1}.
"""
from common import *


def thomas(a, b, c, d):
    """Возвращает (P, Q, x): P[i-1] = P_i, Q[i-1] = Q_i для i = 1..n+1 (всего n+1 значений)."""
    n = len(b)
    B = [-v for v in b]                       # каноническая форма методички
    P, Q = [0.0] * (n + 1), [0.0] * (n + 1)   # P_1 = Q_1 = 0
    for i in range(n):                        # прямой ход: i = 1..n
        denom = B[i] - a[i] * P[i]
        if abs(denom) < EPS:
            raise ValueError(f"Нулевой знаменатель в строке {i + 1}: прогонка неприменима")
        P[i + 1] = c[i] / denom + 0.0          # +0.0 убирает «-0.0»
        Q[i + 1] = (a[i] * Q[i] - d[i]) / denom
    x = [0.0] * n                             # обратный ход
    x[n - 1] = Q[n]                           # x_n = Q_{n+1}
    for i in range(n - 2, -1, -1):
        x[i] = P[i + 1] * x[i + 1] + Q[i + 1]
    return P, Q, x


def main():
    params, _ = read_input(input_path("thomas.txt"))
    a, b, c, d = (parse_vector(params[k]) for k in ("a", "b", "c", "d"))
    n = len(b)
    if not (len(a) == len(c) == len(d) == n):
        raise ValueError("Векторы a, b, c, d должны иметь одинаковую длину n (a[1] и c[n] = 0)")
    a[0], c[n - 1] = 0.0, 0.0

    header("1.2 Метод прогонки")
    print("Нижняя диагональ a:", a)
    print("Главная диагональ b:", b)
    print("Верхняя диагональ c:", c)
    print("Правая часть d:", d)
    dominant = all(abs(b[i]) >= abs(a[i]) + abs(c[i]) for i in range(n))
    print("Диагональное преобладание (достаточное условие устойчивости):", "выполнено" if dominant else "НЕ выполнено")
    print("Каноническая форма методички: a_i x_{i-1} - B_i x_i + c_i x_{i+1} = d_i, где B_i = -b_i:", [-v for v in b])

    P, Q, x = thomas(a, b, c, d)
    print("\n1) Прогоночные коэффициенты (P_1 = Q_1 = 0; x_n = Q_{n+1}):")
    print("   i |          P_i |          Q_i")
    for i in range(n + 1):
        print(f"  {i + 1:2d} | {P[i]:12.8f} | {Q[i]:12.8f}")
    print("\n2) Решение системы:")
    print(fmt_vector(x))
    # проверка по трём диагоналям, без полной матрицы
    res = max(abs((a[i] * x[i - 1] if i > 0 else 0) + b[i] * x[i] + (c[i] * x[i + 1] if i < n - 1 else 0) - d[i])
              for i in range(n))
    print(f"   Невязка max|Ax - d| = {res:.3e}")
    det = 1.0
    for i in range(n):                         # det A = prod(a_i P_i - B_i) для формы методички (с. 23)
        det *= (a[i] * P[i] - (-b[i]))
    print(f"   det A (по прогоночным коэффициентам, методичка с. 23) = {det:.6f}")


if __name__ == "__main__":
    main()
