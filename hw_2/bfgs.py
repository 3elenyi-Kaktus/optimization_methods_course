import sys
from copy import deepcopy

import numpy as np


def f1(x, a, b):
    return (x[0] - x[1] ** 2) ** 2 + (a - x[0]) ** 2 + b


def f2(x, a, b):
    return (x[0] - a) ** 2 + x[0] * x[1] + (x[1] - b) ** 2


def df1_1(x, a, b):
    return -2 * (a - 2 * x[0] + x[1] ** 2)


def df1_2(x, a, b):
    return -4 * x[1] * (x[0] - x[1] ** 2)


def df2_1(x, a, b):
    return -2 * a + 2 * x[0] + x[1]


def df2_2(x, a, b):
    return -2 * b + x[0] + 2 * x[1]


def superf(t_k, coef):
    return f(coef[0] + t_k * coef[1], coef[2], coef[3])


def fib(n):
    if n in [0, 1]:
        return 1
    return fib(n - 1) + fib(n - 2)


def fib_search(f, bounds, tol, coef, max_eps=0.01):
    fib_barrier = (bounds[1] - bounds[0]) / tol
    iter = 0
    while fib(iter) < fib_barrier:
        iter += 1

    k = 1
    left = bounds[0] + fib(iter - k - 1) * (bounds[1] - bounds[0]) / fib(iter - k + 1)
    right = bounds[0] + fib(iter - k) * (bounds[1] - bounds[0]) / fib(iter - k + 1)
    f_left = f(left, coef)
    f_right = f(right, coef)

    while k < iter - 2:
        if f_left <= f_right:
            bounds[1] = right
            right = left
            f_right = f_left
            left = bounds[0] + fib(iter - k - 2) * (bounds[1] - bounds[0]) / fib(iter - k)
            f_left = f(left, coef)
        else:
            bounds[0] = left
            left = right
            f_left = f_right
            right = bounds[0] + fib(iter - k - 1) * (bounds[1] - bounds[0]) / fib(iter - k)
            f_right = f(right, coef)
        k += 1

    right = left + max_eps
    f_right = f(right, coef)
    if f_left >= f_right:
        bounds[0] = left
    else:
        bounds[1] = right
    return (bounds[0] + bounds[1]) / 2


def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)


def ft(t, dk, x, f, a, b):
    x += t * dk
    return f(x, a, b)


def bfgs_method(f, df, a, b, x0, M, t0, eps):
    x_k = [deepcopy(np.array([x0[0], x0[1]]))] * (M + 1)
    for k in range(M):
        grad = [df[0](x_k[k], a, b), df[1](x_k[k], a, b)]
        if np.linalg.norm(grad) < eps:
            return f(x_k[k], a, b)

        if k == 0:
            ck = np.eye(2, dtype=float)

        dk = -ck @ np.array(grad)
        coef = [x_k[k], dk, a, b]
        tk = fib_search(superf, [0, t0], 0.01, coef)

        x_k[k + 1] = x_k[k] + tk * dk

        sk = np.ndarray((2, 1), dtype=float, buffer=x_k[k + 1] - x_k[k])
        yk = np.ndarray((2, 1), dtype=float, buffer=np.array([df[0](x_k[k + 1], a, b), df[1](x_k[k + 1], a, b)]) - np.array(grad))

        pk = 1 / (yk.T @ sk)[0][0]
        ck = (np.eye(2, dtype=float) - (pk * sk) @ yk.T) @ ck @ (np.eye(2, dtype=float) - (pk * yk) @ sk.T) + (
                    pk * sk) @ sk.T
    return f(x_k[M], a, b)


if __name__ == '__main__':
    type = int(input())
    a, b = map(float, input().split(" "))
    if type == 0:
        f = f1
        df = [df1_1, df1_2]
    else:
        f = f2
        df = [df2_1, df2_2]
    eps = float(input())
    x1, x2 = map(float, input().split(" "))
    print(bfgs_method(f, df, a, b, [x1, x2], 10000, 0.5, eps))
