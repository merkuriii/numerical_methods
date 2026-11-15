import os
from common import *


def divided_differences(xs, ys): 
    n = len(xs)
    table = [ys[:]] 
    for k in range(1, n):
        prev = table[-1]
        table.append([(prev[i + 1] - prev[i]) / (xs[i + k] - xs[i]) for i in range(n - k)])
    return table


def lagrange_weights(xs, ys):
    w = []
    for i in range(len(xs)):
        den = 1.0
        for j in range(len(xs)):
            if j != i:
                den *= xs[i] - xs[j]
        w.append(ys[i] / den)
    return w


def lagrange_eval(xs, w, x):
    total = 0.0
    for i in range(len(xs)):
        term = w[i]
        for j in range(len(xs)):
            if j != i:
                term *= x - xs[j]
        total += term
    return total


def lagrange_poly(xs, w):
    result = [0.0]
    for i in range(len(xs)):
        term = [w[i]]
        for j in range(len(xs)):
            if j != i:
                term = poly_mul(term, [-xs[j], 1.0])
        result = poly_add(result, term)
    return result


def newton_coeffs(xs, ys):
    return [row[0] for row in divided_differences(xs, ys)]


def newton_eval(xs, coef, x):
    r = coef[-1]
    for k in range(len(coef) - 2, -1, -1):
        r = r * (x - xs[k]) + coef[k]
    return r


def newton_poly(xs, coef):
    result = [0.0]
    basis = [1.0]
    for k, c in enumerate(coef):
        result = poly_add(result, [c * v for v in basis])
        basis = poly_mul(basis, [-xs[k], 1.0])
    return result


def factor_str(x0):
    return "x" if abs(x0) < 1e-12 else f"(x - {x0:.6f})"


def newton_str(xs, coef):
    s = ""
    for k, c in enumerate(coef):
        factors = "*".join(factor_str(xs[j]) for j in range(k))
        term = f"{abs(c):.6f}" + (f"*{factors}" if factors else "")
        s += (("-" if c < 0 else "") + term) if not s else (" - " if c < 0 else " + ") + term
    return s


def plot(name, f, xs, ys, xstar, L, N, title):
    plt = get_plt()
    if plt is None:
        print("(matplotlib не установлен, график пропущен)")
        return
    lo, hi = min(xs + [xstar]), max(xs + [xstar])
    pad = 0.15 * (hi - lo)
    grid = [lo - pad + (hi - lo + 2 * pad) * i / 400 for i in range(401)]
    plt.figure(figsize=(7.5, 4.8))
    plt.plot(grid, [f(t) for t in grid], color="gray", lw=2.5, label="f(x)")
    plt.plot(grid, [L(t) for t in grid], color="tab:blue", lw=1.5, label="многочлен Лагранжа")
    plt.plot(grid, [N(t) for t in grid], color="tab:orange", lw=1.5, ls="--", label="многочлен Ньютона")
    plt.plot(xs, ys, "ko", label="узлы")
    plt.axvline(xstar, color="tab:red", ls=":", label=f"x* = {xstar:g}")
    plt.plot(xstar, L(xstar), "r*", ms=10)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.legend()
    plt.grid(alpha=0.3)
    save_plot(plt, name)


def main():
    path = input_path("interp_a.txt")
    p = read_params(path)
    consts = constants(p)
    f = compile_expr(p["f"], ["x"], consts)
    xs = expr_list(p["nodes"], consts)
    xstar = float(p["xstar"])
    ys = [f(x) for x in xs]
    tag = os.path.splitext(os.path.basename(path))[0]

    w = lagrange_weights(xs, ys)
    table = divided_differences(xs, ys)
    coef = [row[0] for row in table]
    L, N = lagrange_eval(xs, w, xstar), newton_eval(xs, coef, xstar)
    exact = f(xstar)

    header(f"3.1 Интерполяция ({tag}), f(x) = {p['f']}")
    print("  i         x_i         y_i")
    for i in range(len(xs)):
        print(f"{i:3d}  {xs[i]:10.6f}  {ys[i]:10.6f}")
    print("\nМногочлен Лагранжа: L(x) = sum w_i * prod_{j!=i}(x - x_j), w_i = y_i / prod_{j!=i}(x_i - x_j)")
    print("  w = [" + ", ".join(f"{v + 0.0:.6f}" for v in w) + "]")
    print("  L(x) = " + poly_str(lagrange_poly(xs, w)))
    print("\nТаблица разделённых разностей:")
    print("  i         x_i  " + "  ".join(f"{'f[' + str(k + 1) + ']':>11}" for k in range(len(xs))))
    for i in range(len(xs)):
        row = [table[k][i] for k in range(len(xs) - i)]
        print(f"{i:3d}  {xs[i]:10.6f}  " + "  ".join(f"{v:11.6f}" for v in row))
    print("\nМногочлен Ньютона: N(x) = " + newton_str(xs, coef))
    print("  N(x) = " + poly_str(newton_poly(xs, coef)))
    print(f"\nx* = {xstar}")
    print(f"L(x*) = {L:.10f}")
    print(f"N(x*) = {N:.10f}")
    print(f"f(x*) = {exact:.10f}")
    print(f"абсолютная погрешность |f(x*) - L(x*)| = {abs(exact - L):.3e}")
    print(f"разность |L(x*) - N(x*)| = {abs(L - N):.3e}")
    plot(f"3_1_{tag}.png", f, xs, ys, xstar, lambda t: lagrange_eval(xs, w, t), lambda t: newton_eval(xs, coef, t),
         f"Интерполяция, {tag}: узлы " + ", ".join(f"{x:.3f}" for x in xs))


if __name__ == "__main__":
    main()
