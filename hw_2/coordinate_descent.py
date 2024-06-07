import math
from copy import deepcopy

import numpy as np


def f1(x, a, b):
    return a * math.sin(x[0]) + b * math.cos(x[1])


def f2(x, a, b):
    return (x[0] - a) ** 2 + x[0] * x[1] + (x[1] - b) ** 2


def df1_1(x, a, b):
    return a * math.cos(x[0])


def df1_2(x, a, b):
    return -b * math.sin(x[1])


def df2_1(x, a, b):
    return -2 * a + 2 * x[0] + x[1]


def df2_2(x, a, b):
    return -2 * b + x[0] + 2 * x[1]


def fun(x, a, b):
    return 2*x[0]**2 + x[0] * x[1] + x[1]**2

def dfun1(x, a, b):
    return 4 * x[0] + x[1]

def dfun2(x, a, b):
    return x[0] + 2 * x[1]

def coordinate_descent(f, df, a, b, x00, M, t0, eps1, eps2):
    x_curr = x00
    flags = [False] * len(x00)
    for j in range(M):
        x = [deepcopy(x_curr)] * (len(x00) + 1)
        tk = t0
        for k in range(len(x00)):
            # print("K:", k)
            grad = [df[0](x[k], a, b), df[1](x[k], a, b)]
            if np.linalg.norm(grad) < eps1:
                return f(x[k], a, b)
            basis = np.zeros(len(x00))
            basis[k] = 1.0
            # print("x_jk", x[k])
            # print("df[x_jk]", df[k](x[k], a, b))
            x[k + 1] = x[k] - tk * df[k](x[k], a, b) * basis
            # print(x[k + 1])
            # print()
            if f(x[k + 1], a, b) - f(x[k], a, b) >= 0:
                tk /= 2
                flags[k] = False
            else:
                if np.linalg.norm(x[k + 1] - x[k]) < eps2 and np.abs(f(x[k + 1], a, b) - f(x[k], a, b)) < eps2:
                    if flags[k]:
                        return f(x[k + 1], a, b)
                    else:
                        flags[k] = True
                else:
                    flags[k] = False
        x_curr = x[len(x00)]
    return f(x_curr, a, b)


if __name__ == "__main__":
    type = int(input())
    a, b = map(float, input().split(" "))
    if type == 0:
        f = f1
        df = [df1_1, df1_2]
    else:
        f = f2
        df = [df2_1, df2_2]
    eps1, eps2 = map(float, input().split(" "))
    x1, x2 = map(float, input().split(" "))
    x = [x1, x2]
    # f = fun
    # df = [dfun1, dfun2]
    # a = 0
    # b = 0
    print(coordinate_descent(f, df, a, b, x, 10000, 0.5, eps1, eps2))
