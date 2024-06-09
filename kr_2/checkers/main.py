import numpy as np
import scipy.optimize as scopt


def main():
    c = -np.array(list(map(float, input().split(" "))))
    a = list(map(float, input().split(" ")))
    b = list(map(float, input().split(" ")))
    b_eq = np.array([a[-1], b[-1]])
    a.pop()
    b.pop()
    A_eq = np.array([a,
                     b,])
    # проставляем индикаторы для типа переменных, 0 - непрерывное, 1 - целое число
    var_types = [2, 2, 2, 2]
    # также указываем границы, в том числе и для целочисленных переменных
    bounds = [(0, None), (0, None), (0, None), (0, None)]
    res_milp = scopt.linprog(c=c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, integrality=var_types)
    print(f"Solution: x = {list(np.round(res_milp['x'], 3))}")
    print(f"f = {-res_milp['fun']}")
    print(f"{res_milp['message']}")


if __name__ == '__main__':
    main()


# -3 -1 -5 4
# 4 0 1 4 2
# -4 1 0 4 1