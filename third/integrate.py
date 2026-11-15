from common import *

ORDERS = {"левые прямоугольники": 1, "правые прямоугольники": 1, "средние прямоугольники": 2,
          "трапеции": 2, "Симпсон": 4}


def steps(a, b, h):
    n = round((b - a) / h)
    if abs(n * h - (b - a)) > 1e-9:
        raise ValueError(f"Шаг {h} не укладывается целое число раз в [{a}, {b}]")
    return n


def left(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + i * h) for i in range(n))


def right(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + i * h) for i in range(1, n + 1))


def middle(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + (i + 0.5) * h) for i in range(n))


def trapezoid(f, a, b, n):
    h = (b - a) / n
    return h * ((f(a) + f(b)) / 2 + sum(f(a + i * h) for i in range(1, n)))


def simpson(f, a, b, n):
    if n % 2:
        raise ValueError("Для метода Симпсона нужно чётное число отрезков")
    h = (b - a) / n
    s = f(a) + f(b)
    s += 4 * sum(f(a + i * h) for i in range(1, n, 2))
    s += 2 * sum(f(a + i * h) for i in range(2, n, 2))
    return h * s / 3


METHODS = {"левые прямоугольники": left, "правые прямоугольники": right, "средние прямоугольники": middle,
           "трапеции": trapezoid, "Симпсон": simpson}


def runge_romberg(i1, i2, h1, h2, p): # h1 > h2
    corr = (i2 - i1) / ((h1 / h2) ** p - 1)
    return i2 + corr, abs(corr) 


def main():
    p = read_params(input_path("integral.txt"))
    consts = constants(p)
    f = compile_expr(p["f"], ["x"], consts)
    a, b, h1, h2 = (float(p[k]) for k in ("a", "b", "h1", "h2"))
    n1, n2 = steps(a, b, h1), steps(a, b, h2)
    exact = None
    if "F" in p:
        F = compile_expr(p["F"], ["x"], consts)
        exact = F(b) - F(a)

    header(f"3.5 Определённый интеграл от {p['f']} на [{a}, {b}]")
    print(f"шаги h1 = {h1} (n = {n1}), h2 = {h2} (n = {n2})")
    cols = f"{'метод':24}{'I(h1)':>12}{'I(h2)':>12}{'p':>3}{'I Рунге-Ромберг':>17}{'оценка погр.(h2)':>18}"
    if exact is not None:
        cols += f"{'|I(h2) - точное|':>18}"
    print(cols)
    rows = []
    for name, method in METHODS.items():
        i1, i2 = method(f, a, b, n1), method(f, a, b, n2)
        order = ORDERS[name]
        refined, est = runge_romberg(i1, i2, h1, h2, order)
        line = f"{name:24}{i1:12.6f}{i2:12.6f}{order:3d}{refined:17.6f}{est:18.2e}"
        err = None
        if exact is not None:
            err = abs(i2 - exact)
            line += f"{err:18.2e}"
        print(line)
        rows.append((name, err, abs(refined - exact) if exact is not None else None))
    if exact is not None:
        print(f"\nточное значение интеграла = {exact:.8f}")
        best = min(rows, key=lambda r: r[1])
        print(f"наиболее точный метод (шаг h2): {best[0]}, погрешность {best[1]:.2e}")
        print("погрешность после Рунге-Ромберга: " + "; ".join(f"{n}: {e:.1e}" for n, _, e in rows))


if __name__ == "__main__":
    main()
