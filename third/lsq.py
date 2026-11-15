from common import *


def fit(xs, ys, m):
    A = [[sum(x ** (k + i) for x in xs) for i in range(m + 1)] for k in range(m + 1)]
    b = [sum(y * x ** k for x, y in zip(xs, ys)) for k in range(m + 1)]
    return solve(A, b)


def errors(xs, ys, coef):
    sse = sum((poly_eval(coef, x) - y) ** 2 for x, y in zip(xs, ys))
    return sse, sse / len(xs)


def plot(xs, ys, fits):
    plt = get_plt()
    if plt is None:
        print("(matplotlib не установлен, график пропущен)")
        return
    lo, hi = min(xs), max(xs)
    pad = 0.1 * (hi - lo)
    grid = [lo - pad + (hi - lo + 2 * pad) * i / 300 for i in range(301)]
    plt.figure(figsize=(7.5, 4.8))
    plt.plot(xs, ys, "ko", label="исходные данные")
    for m, coef in fits:
        plt.plot(grid, [poly_eval(coef, t) for t in grid], lw=1.8, label=f"многочлен степени {m}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Метод наименьших квадратов")
    plt.legend()
    plt.grid(alpha=0.3)
    save_plot(plt, "3_3_lsq.png")


def main():
    p = read_params(input_path("lsq.txt"))
    xs, ys = vec(p["x"]), vec(p["y"])
    degrees = [int(t) for t in p["degrees"].split()]
    fits = [(m, fit(xs, ys, m)) for m in degrees]

    header("3.3 Метод наименьших квадратов")
    print("  i         x_i         y_i")
    for i in range(len(xs)):
        print(f"{i:3d}  {xs[i]:10.6f}  {ys[i]:10.6f}")
    res = []
    for m, coef in fits:
        sse, mse = errors(xs, ys, coef)
        res.append((m, sse, mse))
        print(f"\nстепень {m}: коэффициенты a_0..a_{m} = [" + ", ".join(f"{v:.6f}" for v in coef) + "]")
        print(f"  P{m}(x) = " + poly_str(coef))
        print(f"  сумма квадратов ошибок = {sse:.6f},  MSE{m} = {mse:.6f}")
    print("\nсравнение: " + ", ".join(f"MSE{m} = {mse:.6f}" for m, _, mse in res))
    best = min(res, key=lambda t: t[2])
    print(f"наименьшая ошибка у многочлена степени {best[0]} (MSE = {best[2]:.6f})")
    plot(xs, ys, fits)


if __name__ == "__main__":
    main()
