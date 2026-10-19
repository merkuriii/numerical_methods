import math
from common import *


def off_norm(A):
    n = len(A)
    return math.sqrt(sum(A[i][j] ** 2 for i in range(n) for j in range(n) if i != j))


def jacobi(A, eps, max_iter=100000):
    n = len(A)
    A = [r[:] for r in A]
    V = identity(n)
    history = [(0, None, None, off_norm(A))]
    for k in range(1, max_iter + 1):
        if off_norm(A) < eps:
            return [A[t][t] for t in range(n)], V, k - 1, history
        i, j = max(((i, j) for i in range(n) for j in range(i + 1, n)), key=lambda p: abs(A[p[0]][p[1]]))
        phi = 0.5 * math.atan2(2 * A[i][j], A[i][i] - A[j][j])
        c, s = math.cos(phi), math.sin(phi)
        U = identity(n)
        U[i][i] = U[j][j] = c
        U[i][j], U[j][i] = -s, s
        A = matmul(matmul(transpose(U), A), U)
        V = matmul(V, U)
        history.append((k, i + 1, j + 1, off_norm(A)))
    raise RuntimeError("Метод не сошёлся")


def main():
    params, A = read_input(input_path("jacobi.txt"))
    n = len(A)
    if any(len(r) != n for r in A):
        raise ValueError("Матрица должна быть квадратной")
    if any(abs(A[i][j] - A[j][i]) > 1e-12 for i in range(n) for j in range(n)):
        raise ValueError("Нужна симметричная матрица")
        
    eps = float(params.get("eps", "1e-8"))
    vals, V, iters, _ = jacobi(A, eps)
    Lam = [[vals[i] if i == j else 0.0 for j in range(n)] for i in range(n)]

    print("1.4 Метод вращений Якоби")
    print(f"eps = {eps}")
    print(fmt_vector(vals, "\nСобственные значения:"))
    print(fmt_matrix(V, "\nМатрица собственных векторов V:", 6))
    
    AV, VL = matmul(A, V), matmul(V, Lam)
    print(fmt_matrix(AV, "\nA * V:", 6))
    print(fmt_matrix(VL, "V * Lambda:", 6))
    print(f"max|A*V - V*Lambda| = {max_abs(mat_sub(AV, VL)):.2e}")


if __name__ == "__main__":
    main()
