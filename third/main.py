import os
import sys

import interpolation, spline, lsq, derivative, integrate

INPUTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs")

TASKS = {
    "3.1": [(interpolation, "interp_a.txt"), (interpolation, "interp_b.txt")],
    "3.2": [(spline, "spline.txt")],
    "3.3": [(lsq, "lsq.txt")],
    "3.4": [(derivative, "derivative.txt")],
    "3.5": [(integrate, "integral.txt")],
}


def run(module, name):
    sys.argv = [module.__name__ + ".py", os.path.join(INPUTS, name)]
    module.main()


def main():
    task = sys.argv[1] if len(sys.argv) > 1 else "all"
    if task != "all" and task not in TASKS:
        print("python main.py [3.1|3.2|3.3|3.4|3.5]")
        return
    for t, jobs in TASKS.items():
        if task in ("all", t):
            for module, name in jobs:
                run(module, name)


if __name__ == "__main__":
    main()
