import math
from copy import deepcopy

import numpy as np


def f1(x, a, b):
    return (x[0] - a) ** 2 + x[1] ** 2 + x[0] / (math.fabs(x[1]) + b)


def f2(x, a, b):
    return (x[0] - a) ** 2 + x[0] * x[1] + (x[1] - b) ** 2


def df1_1(x, a, b):
    return -2 * a + 1 / (math.fabs(x[1]) + b) + 2 * x[0]


def df1_2(x, a, b):
    return 2 * x[1] - (x[0] * x[1]) / (math.fabs(x[1]) * (math.fabs(x[1]) + b) ** 2)


def df2_1(x, a, b):
    return -2 * a + 2 * x[0] + x[1]


def df2_2(x, a, b):
    return -2 * b + x[0] + 2 * x[1]


def fib_1(t_k, coef):
    return f1(coef[0] + t_k * coef[1], coef[2], coef[3])


def fib_2(t_k, coef):
    return f2(coef[0] + t_k * coef[1], coef[2], coef[3])


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


def fletcher_reeves(f, df, a, b, x0, M, t0, eps1, eps2):
    x_k = [deepcopy(np.array([x0[0], x0[1]]))] * (M + 1)
    flag = False
    for k in range(M):
        grad = [df[0](x_k[k], a, b), df[1](x_k[k], a, b)]
        if np.linalg.norm(grad) < eps1:
            return f(x_k[k], a, b)
        if k == 0:
            dk = np.array([-df[0](x0, a, b), -df[1](x0, a, b)])

            coef = [x_k[k], dk, a, b]
            t_k = fib_search(fibo, [0, 0.5], 0.01, coef)

            x_k[k + 1] = x_k[k] + t_k * dk

            if np.linalg.norm(x_k[k + 1] - x_k[k]) < eps2 and np.abs(f(x_k[k + 1], a, b) - f(x_k[k], a, b)) < eps2:
                if flag:
                    return f(x_k[k + 1], a, b)
                flag = True
            else:
                flag = False
        else:
            grad_curr = [df[0](x_k[k], a, b), df[1](x_k[k], a, b)]
            grad_prev = [df[0](x_k[k - 1], a, b), df[1](x_k[k - 1], a, b)]
            beta = (np.linalg.norm(grad_curr) ** 2) / (np.linalg.norm(grad_prev) ** 2)
            dk = np.array([-grad_curr[0], -grad_curr[1]]) + beta * dk

            coef = [x_k[k], dk, a, b]
            t_k = fib_search(fibo, [0, 0.5], 0.01, coef)
            # print(t_k)

            x_k[k + 1] = x_k[k] + t_k * dk

            if np.linalg.norm(x_k[k + 1] - x_k[k]) < eps2 and np.abs(f(x_k[k + 1], a, b) - f(x_k[k], a, b)) < eps2:
                if flag:
                    return f(x_k[k + 1], a, b)
                flag = True
            else:
                flag = False
    return f(x_k[M - 1], a, b)


if __name__ == '__main__':
    type = int(input())
    a, b = map(float, input().split(" "))
    if type == 0:
        f = f1
        df = [df1_1, df1_2]
        fibo = fib_1
    else:
        f = f2
        df = [df2_1, df2_2]
        fibo = fib_2
    eps1, eps2 = map(float, input().split(" "))
    x1, x2 = map(float, input().split(" "))
    print(fletcher_reeves(f, df, a, b, [x1, x2], 10000, 0.5, eps1, eps2))
