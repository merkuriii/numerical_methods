import os
import sys

import lu, thomas, iterative, jacobi, qr_alg

TASKS = {
    "1.1": (lu, "lu.txt"),
    "1.2": (thomas, "thomas.txt"),
    "1.3": (iterative, "iterative.txt"),
    "1.4": (jacobi, "jacobi.txt"),
    "1.5": (qr_alg, "qr.txt"),
}
INPUTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs")


def run(task, path=None):
    module, default = TASKS[task]
    sys.argv = [module.__name__ + ".py", path or os.path.join(INPUTS, default)]
    module.main()


def main():
    args = sys.argv[1:]
    if not args:
        for t in TASKS:
            run(t)
        print("\n1.5, пример с комплексными собственными значениями")
        run("1.5", os.path.join(INPUTS, "qr_complex.txt"))
    elif args[0] in TASKS:
        run(args[0], args[1] if len(args) > 1 else None)
    else:
        print("python main.py [1.1|1.2|1.3|1.4|1.5] [файл]")


if __name__ == "__main__":
    main()
