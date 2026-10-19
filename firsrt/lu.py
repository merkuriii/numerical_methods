from common import *


def lu_decompose(A):
    n = len(A)
    U = [r[:] for r in A] 
    L = [[0.0] * n for _ in range(n)]
    perm = list(range(n))
    sign = 1
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(U[i][k]))
        if abs(U[p][k]) < EPS:
            raise ValueError("Матрица вырождена: определитель равен 0")
        if p != k:
            U[k], U[p] = U[p], U[k]
            L[k], L[p] = L[p], L[k]
            perm[k], perm[p] = perm[p], perm[k]
            sign = -sign
        for i in range(k + 1, n):
            m = U[i][k] / U[k][k]
            L[i][k] = m
            for j in range(k, n):
                U[i][j] -= m * U[k][j]
    for i in range(n):
        L[i][i] = 1.0
    return L, U, perm, sign


def perm_matrix(perm):
    n = len(perm)
    return [[1.0 if j == perm[i] else 0.0 for j in range(n)] for i in range(n)]


def lu_solve(L, U, perm, b):
    n = len(L)
    y = [0.0] * n
    for i in range(n):                       # L*y = P*b
        y[i] = b[perm[i]] - sum(L[i][j] * y[j] for j in range(i))
    x = [0.0] * n
    for i in range(n - 1, -1, -1):           # U*x = y
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x


def determinant(U, sign):
    d = float(sign)
    for i in range(len(U)):
        d *= U[i][i]
    return d


def inverse(L, U, perm):
    n = len(L)
    cols = []
    for j in range(n):
        e = [1.0 if i == j else 0.0 for i in range(n)]
        cols.append(lu_solve(L, U, perm, e))
    return transpose(cols)


def main():
    _, rows = read_input(input_path("lu.txt"))
    A, b = split_augmented(rows)
    n = len(A)
    L, U, perm, sign = lu_decompose(A)
    P = perm_matrix(perm)
    LU = matmul(L, U)
    PA = matmul(P, A)
    x = lu_solve(L, U, perm, b)
    Ainv = inverse(L, U, perm)
    AAinv = matmul(A, Ainv)
    det_A = determinant(U, sign)

    print(fmt_matrix(L, "L:"))
    print(fmt_matrix(U, "\nU:"))
    print(fmt_matrix(LU, "\nL * U:"))
    print("\nРешение x:")
    print(fmt_vector(x))
    print(fmt_matrix(Ainv, "\nA^-1:"))
    print(f"\nDeterminant det(A) = {det_A:.6f}")
    print(fmt_matrix(AAinv, "\nA * A^-1:"))
    print(f"max|A * A^-1 - E| = {max_abs(mat_sub(AAinv, identity(n))):.2e}")
    print(fmt_matrix(PA, "\nP * A:"))
    print(f"max|L * U - P * A| = {max_abs(mat_sub(LU, PA)):.2e}")


if __name__ == "__main__":
    main()
