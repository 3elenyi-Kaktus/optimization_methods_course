import numpy as np


def f0(x, coeff):
    return coeff[0] * x[0] ** 4 + coeff[1] * x[1] ** 3 + coeff[2] * x[1] ** 2 + coeff[3] * x[0] + coeff[4]


def f1(x, coeff):
    return x[0] ** 2 + coeff[0] * x[0] * x[1] + coeff[1] * (x[1] - 3) ** 2


def jeeves_hooke(f, x0, tolerance, coeff):
    delta = 1.0
    al = 2.0
    x_best = x0.copy()
    x_curr = x0.copy()
    while delta >= tolerance:
        for pos, coord in enumerate(x_curr):
            x = x_curr.copy()
            x[pos] += delta
            if f(x, coeff) < f(x_best, coeff):
                x_curr = x
                continue
            x[pos] -= 2 * delta
            if f(x, coeff) < f(x_best, coeff):
                x_curr = x
        if f(x_curr, coeff) >= f(x_best, coeff):
            delta /= al
            x_curr = x_best.copy()
        else:
            x_best = x_curr.copy()
            x_curr = 2 * x_curr - x_best
    return x_best


type_ = int(input())
f_ = f0 if (type_ == 0) else f1
coef = [i for i in map(float, input().split())]
x0_ = np.array([i for i in map(float, input().split())])
tol = float(input())
r1 = jeeves_hooke(f_, x0_, tol, coef)
print("{:.10f} {:.10f}".format(r1[0], r1[1]))
