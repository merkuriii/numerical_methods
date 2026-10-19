import cmath
import math
from common import *


def qr_decompose(A):
    n = len(A)
    R = [r[:] for r in A]
    Q = identity(n)
    for k in range(n - 1):
        x = [R[i][k] for i in range(k, n)]
        norm = math.sqrt(sum(t * t for t in x))
        if norm < EPS:
            continue
        alpha = -norm if x[0] >= 0 else norm
        v = x[:]
        v[0] -= alpha
        vn = math.sqrt(sum(t * t for t in v))
        if vn < EPS:
            continue
        v = [t / vn for t in v]
        H = identity(n)
        for i in range(len(v)):
            for j in range(len(v)):
                H[k + i][k + j] -= 2 * v[i] * v[j]
        R = matmul(H, R)
        Q = matmul(Q, H)
    for i in range(n):
        for j in range(i):
            R[i][j] = 0.0
    return Q, R


def eigenvalues_from_form(A, eps):
    n, vals, i = len(A), [], 0
    while i < n:
        if i == n - 1 or abs(A[i + 1][i]) <= eps:
            vals.append(complex(A[i][i]))
            i += 1
        else:
            a, b, c, d = A[i][i], A[i][i + 1], A[i + 1][i], A[i + 1][i + 1]
            tr, det = a + d, a * d - b * c
            disc = cmath.sqrt(tr * tr - 4 * det)
            vals += [(tr + disc) / 2, (tr - disc) / 2]
            i += 2
    return vals


def qr_algorithm(A, eps, max_iter=100000):
    n = len(A)
    Ak = [r[:] for r in A]
    prev = None
    for k in range(1, max_iter + 1):
        Q, R = qr_decompose(Ak)
        Ak = matmul(R, Q)
        below = max([abs(Ak[i][j]) for i in range(n) for j in range(i - 1)] + [0.0])
        vals = sorted(eigenvalues_from_form(Ak, eps), key=lambda z: (z.real, z.imag))
        if prev is not None and below < eps and max(abs(a - b) for a, b in zip(vals, prev)) < eps:
            return vals, Ak, k
        prev = vals
    raise RuntimeError("QR-алгоритм не сошёлся")


def fmt_c(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:.8f}"
    return f"{z.real:.8f} {'+' if z.imag >= 0 else '-'} {abs(z.imag):.8f}i"


def main():
    params, A = read_input(input_path("qr.txt"))
    n = len(A)
    if any(len(r) != n for r in A):
        raise ValueError("Матрица должна быть квадратной")
    eps = float(params.get("eps", "1e-9"))
    
    Q, R = qr_decompose(A)
    QR = matmul(Q, R)
    vals, _, _ = qr_algorithm(A, eps)

    print("1.5 QR-разложение")
    print(fmt_matrix(A, "1. Исходная матрица A:"))
    print(fmt_matrix(Q, "\n2. Матрица Q:"))
    print(fmt_matrix(R, "\nМатрица R:"))
    print(fmt_matrix(QR, "\n3. Результат умножения Q · R:"))
    
    print("\n4. Найденные собственные значения:")
    for i, v in enumerate(vals):
        print(f"  lambda{i + 1} = {fmt_c(v)}")


if __name__ == "__main__":
    main()
