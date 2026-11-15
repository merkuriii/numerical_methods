from common import *


def build_spline(xs, ys):
    n = len(xs) - 1 # отрезков = n, узлов = n + 1
    h = [xs[i + 1] - xs[i] for i in range(n)] # шаги между узлами
    c = [0.0] * (n + 1)
    if n > 1:
        A = [h[j - 1] for j in range(1, n)]
        B = [2 * (h[j - 1] + h[j]) for j in range(1, n)]
        C = [h[j] for j in range(1, n)]
        D = [3 * ((ys[j + 1] - ys[j]) / h[j] - (ys[j] - ys[j - 1]) / h[j - 1]) for j in range(1, n)]
        c[1:n] = tridiagonal(A, B, C, D)
    a = ys[:n]
    b = [(ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3 for j in range(n)]
    d = [(c[j + 1] - c[j]) / (3 * h[j]) for j in range(n)]
    return a, b, c[:n], d, h


def find_interval(xs, x):
    for j in range(len(xs) - 1):
        if xs[j] <= x <= xs[j + 1]:
            return j
    raise ValueError("Точка вне отрезка интерполяции")


def spline_eval(xs, coefs, x):
    a, b, c, d, _ = coefs
    j = find_interval(xs, x)
    t = x - xs[j]
    return a[j] + b[j] * t + c[j] * t ** 2 + d[j] * t ** 3


def second_derivative(xs, coefs, x):
    _, _, c, d, _ = coefs
    j = find_interval(xs, x)
    return 2 * c[j] + 6 * d[j] * (x - xs[j])


def plot(xs, ys, coefs, xstar, sx):
    plt = get_plt()
    if plt is None:
        print("(matplotlib не установлен, график пропущен)")
        return
    plt.figure(figsize=(7.5, 4.8))
    for j in range(len(xs) - 1):
        grid = [xs[j] + (xs[j + 1] - xs[j]) * i / 100 for i in range(101)]
        plt.plot(grid, [spline_eval(xs, coefs, t) for t in grid], lw=2, label=f"участок {j + 1}")
    plt.plot(xs, ys, "ko", label="узлы")
    plt.plot(xstar, sx, "r*", ms=12, label=f"(x*, S(x*)) = ({xstar:g}; {sx:.4f})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Кубический сплайн")
    plt.legend()
    plt.grid(alpha=0.3)
    save_plot(plt, "3_2_spline.png")


def main():
    p = read_params(input_path("spline.txt"))
    xs, ys, xstar = vec(p["x"]), vec(p["y"]), float(p["xstar"])
    coefs = build_spline(xs, ys)
    a, b, c, d, h = coefs
    n = len(xs) - 1
    sx = spline_eval(xs, coefs, xstar)
    k = find_interval(xs, xstar)

    header("3.2 Кубический сплайн (нулевая кривизна на концах)")
    print("  i         x_i         y_i")
    for i in range(n + 1):
        print(f"{i:3d}  {xs[i]:10.6f}  {ys[i]:10.6f}")
    print("\nS_i(x) = a_i + b_i(x - x_{i-1}) + c_i(x - x_{i-1})^2 + d_i(x - x_{i-1})^3 на [x_{i-1}, x_i]")
    print("  i  [x_{i-1}, x_i]            a_i         b_i         c_i         d_i")
    for j in range(n):
        print(f"{j + 1:3d}  [{xs[j]:.4f}, {xs[j + 1]:.4f}]  {a[j]:11.6f} {b[j]:11.6f} {c[j]:11.6f} {d[j]:11.6f}")
    print(f"\nx* = {xstar}: интервал номер {k + 1} = [{xs[k]}, {xs[k + 1]}]")
    print(f"S(x*) = {sx:.8f}")
    dev = max(abs(spline_eval(xs, coefs, xs[i]) - ys[i]) for i in range(n + 1))
    print(f"проверка S(x_i) = y_i: max отклонение = {dev:.2e}")
    print(f"S''(x_0) = {second_derivative(xs, coefs, xs[0]):.2e},  S''(x_n) = {second_derivative(xs, coefs, xs[-1]):.2e}")
    plot(xs, ys, coefs, xstar, sx)


if __name__ == "__main__":
    main()
