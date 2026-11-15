from common import *


def choose_nodes(xs, xstar):
    best = None
    for j in range(len(xs) - 2):
        if xs[j] <= xstar <= xs[j + 2]:
            dist = abs(xstar - xs[j + 1])
            if best is None or dist < best[0] - 1e-12:
                best = (dist, j)
    if best is None:
        raise ValueError("Точка вне таблицы или узлов меньше трёх")
    return best[1]


def derivatives(xs, ys, xstar):
    j = choose_nodes(xs, xstar)
    x0, x1, x2 = xs[j:j + 3]
    y0, y1, y2 = ys[j:j + 3]
    f01 = (y1 - y0) / (x1 - x0)
    f12 = (y2 - y1) / (x2 - x1)
    f012 = (f12 - f01) / (x2 - x0)
    d1 = f01 + f012 * (2 * xstar - x0 - x1)
    d2 = 2 * f012
    return j, d1, d2


def main():
    p = read_params(input_path("derivative.txt"))
    xs, ys, xstar = vec(p["x"]), vec(p["y"]), float(p["xstar"])
    j, d1, d2 = derivatives(xs, ys, xstar)

    header("3.4 Численное дифференцирование")
    print("  i         x_i         y_i")
    for i in range(len(xs)):
        print(f"{i:3d}  {xs[i]:10.6f}  {ys[i]:10.6f}")
    print(f"\nточка x* = {xstar}")
    print("узлы, использованные при вычислении (многочлен Ньютона 2-й степени): " +
          ", ".join(f"x_{j + k} = {xs[j + k]}" for k in range(3)))
    print(f"f'(x*)  = {d1:.6f}")
    print(f"f''(x*) = {d2:.6f}")


if __name__ == "__main__":
    main()
