from common import *
from lu import lu_decompose, lu_solve

MAX_ITER = 100000


def to_iteration_form(A, b):
    n = len(A)
    if any(abs(A[i][i]) < EPS for i in range(n)):
        raise ValueError("Ноль на диагонали, переставьте уравнения")
    alpha = [[0.0 if i == j else -A[i][j] / A[i][i] for j in range(n)] for i in range(n)]
    beta = [b[i] / A[i][i] for i in range(n)]
    return alpha, beta


def simple_iteration(alpha, beta, eps, exact):
    n = len(beta)
    q = mat_norm_inf(alpha)
    factor = q / (1 - q) if q < 1 else 1.0
    x, hist = beta[:], []
    for k in range(1, MAX_ITER + 1):
        x_new = [sum(alpha[i][j] * x[j] for j in range(n)) + beta[i] for i in range(n)]
        err = factor * vec_norm_inf(vec_sub(x_new, x))
        x = x_new
        hist.append(vec_norm_inf(vec_sub(x, exact)))
        if err <= eps:
            return x, k, hist, True
    return x, MAX_ITER, hist, False


def seidel(alpha, beta, eps, exact):
    n = len(beta)
    q = mat_norm_inf(alpha)
    C = [[alpha[i][j] if j > i else 0.0 for j in range(n)] for i in range(n)]
    factor = mat_norm_inf(C) / (1 - q) if q < 1 else 1.0
    x, hist = beta[:], []
    for k in range(1, MAX_ITER + 1):
        x_new = x[:]
        for i in range(n):
            x_new[i] = sum(alpha[i][j] * x_new[j] for j in range(n)) + beta[i]
        err = factor * vec_norm_inf(vec_sub(x_new, x))
        x = x_new
        hist.append(vec_norm_inf(vec_sub(x, exact)))
        if err <= eps:
            return x, k, hist, True
    return x, MAX_ITER, hist, False


def main():
    params, rows = read_input(input_path("iterative.txt"))
    A, b = split_augmented(rows)
    eps = float(params.get("eps", "1e-6"))
    L, U, perm, _ = lu_decompose(A)
    exact = lu_solve(L, U, perm, b)
    alpha, beta = to_iteration_form(A, b)
    _, ks, _, _ = simple_iteration(alpha, beta, eps, exact)
    _, kz, _, _ = seidel(alpha, beta, eps, exact)

    print(f"eps = {eps}")
    print(fmt_vector(exact, "Эталонное решение:"))
    print(f"Количество итераций (простые): {ks}")
    print(f"Количество итераций (Зейдель): {kz}")
    print(f"Сравнение скорости: Зейдель быстрее в {ks / kz:.2f} раза")


if __name__ == "__main__":
    main()
