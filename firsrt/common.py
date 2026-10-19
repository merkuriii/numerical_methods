import os
import sys

EPS = 1e-12


def input_path(name):
    if len(sys.argv) > 1:
        return sys.argv[1]
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs", name)


def parse_vector(s):
    return [float(x) for x in s.replace(",", " ").replace("|", " ").split()]


def read_input(path):
    params, rows = {}, []
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        if "=" in line:
            k, v = line.split("=", 1)
            params[k.strip()] = v.strip()
        else:
            rows.append(parse_vector(line))
    return params, rows


def split_augmented(rows):
    n = len(rows)
    if any(len(r) != n + 1 for r in rows):
        raise ValueError("Нужна матрица n x (n+1): коэффициенты и правая часть в каждой строке")
    return [r[:-1] for r in rows], [r[-1] for r in rows]


def identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(r) for r in zip(*A)]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def mat_sub(A, B):
    return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def vec_norm_inf(v):
    return max(abs(x) for x in v)


def mat_norm_inf(A):
    return max(sum(abs(x) for x in row) for row in A)


def max_abs(A):
    return max(abs(x) for row in A for x in row)


def fmt_matrix(M, title=None, prec=4):
    w = max(len(f"{x:.{prec}f}") for row in M for x in row)
    lines = [title] if title else []
    lines += ["  " + "  ".join(f"{x:{w}.{prec}f}" for x in row) for row in M]
    return "\n".join(lines)


def fmt_vector(v, name="x", prec=6):
    return f"{name} = [" + ", ".join(f"{x:.{prec}f}" for x in v) + "]"
