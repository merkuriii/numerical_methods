import math
import os
import sys

EPS = 1e-12

NAMES = {k: v for k, v in vars(math).items() if not k.startswith("_")}
NAMES.update(ln=math.log, lg=math.log10, tg=math.tan, abs=abs)

HERE = os.path.dirname(os.path.abspath(__file__))


def input_path(name):
    if len(sys.argv) > 1:
        return sys.argv[1]
    return os.path.join(HERE, "inputs", name)


def read_params(path):
    params = {}
    for line in open(path, encoding="utf-8"):
        if "=" in line:
            k, v = line.split("=", 1)
            params[k.strip()] = v.strip()
    return params


def vec(s):
    return [float(t) for t in s.replace(",", " ").split()]


def constants(params):
    out = {}
    for k, v in params.items():
        try:
            out[k] = float(v)
        except ValueError:
            pass
    return out


def compile_expr(expr, args, consts):
    code = compile(expr.strip().replace("^", "**"), "<expr>", "eval")
    env = dict(NAMES)
    env.update(consts)

    def fn(*vals):
        env.update(zip(args, vals))
        return eval(code, {"__builtins__": {}}, env)
    return fn


def expr_list(s, consts):
    return [compile_expr(t, [], consts)() for t in s.split(";")]


def solve(A, b):
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[p][k]) < EPS:
            raise ValueError("Вырожденная система")
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= m * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


def tridiagonal(a, b, c, d):
    n = len(b)
    P, Q = [0.0] * n, [0.0] * n
    for i in range(n):
        den = b[i] + (a[i] * P[i - 1] if i else 0.0)
        P[i] = -c[i] / den if i < n - 1 else 0.0
        Q[i] = (d[i] - (a[i] * Q[i - 1] if i else 0.0)) / den
    x = [0.0] * n
    x[-1] = Q[-1]
    for i in range(n - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]
    return x


def poly_mul(p, q):
    r = [0.0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def poly_add(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0.0) + (q[i] if i < len(q) else 0.0) for i in range(n)]


def poly_eval(c, x):
    r = 0.0
    for a in reversed(c):
        r = r * x + a
    return r


def poly_str(c, prec=6):
    s = ""
    for k, v in enumerate(c):
        if abs(v) < 1e-12:
            continue
        mag = f"{abs(v):.{prec}f}"
        var = "" if k == 0 else ("x" if k == 1 else f"x^{k}")
        term = mag if k == 0 else f"{mag}*{var}"
        s += (("-" if v < 0 else "") + term) if not s else (" - " if v < 0 else " + ") + term
    return s or "0"


def header(text):
    print(f"\n{text}")


def get_plt():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None


def save_plot(plt, name):
    os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)
    plt.savefig(os.path.join(HERE, "plots", name), dpi=130, bbox_inches="tight")
    plt.close()
    print(f"график сохранён: plots/{name}")
