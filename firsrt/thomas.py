from common import *


def thomas(a, b, c, d):
    n = len(b)
    B = [-v for v in b]
    P, Q = [0.0] * (n + 1), [0.0] * (n + 1)
    for i in range(n):                        # прямой ход
        denom = B[i] - a[i] * P[i]
        if abs(denom) < EPS:
            raise ValueError(f"Нулевой знаменатель в строке {i + 1}: прогонка неприменима")
        P[i + 1] = c[i] / denom + 0.0
        Q[i + 1] = (a[i] * Q[i] - d[i]) / denom
    x = [0.0] * n                             # обратный ход
    x[n - 1] = Q[n]
    for i in range(n - 2, -1, -1):
        x[i] = P[i + 1] * x[i + 1] + Q[i + 1]
    return P, Q, x


def main():
    params, _ = read_input(input_path("thomas.txt"))
    a, b, c, d = (parse_vector(params[k]) for k in ("a", "b", "c", "d"))
    n = len(b)
    if not (len(a) == len(c) == len(d) == n):
        raise ValueError("Векторы a, b, c, d должны иметь одинаковую длину n")
    a[0], c[n - 1] = 0.0, 0.0

    P, Q, x = thomas(a, b, c, d)

    print("1) Прогоночные коэффициенты:")
    print("   i |          P_i |          Q_i")
    for i in range(n + 1):
        print(f"   {i + 1:2d} | {P[i]:12.8f} | {Q[i]:12.8f}")

    print("\n2) Решение системы:")
    print(fmt_vector(x))


if __name__ == "__main__":
    main()
