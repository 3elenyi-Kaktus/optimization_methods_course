import numpy as np


def F(x, coef, r):
    return f(x, coef) + r / 2 * (g(x, coef) ** 2 + (x[0] < 0) * (x[0] ** 2) + (x[1] < 0) * (x[1] ** 2))


def f(x, coef):
    return -coef[0] * x[0] + coef[1] * (x[1] ** 2)


def g(x, coef):
    return x[0] ** 3 - x[1] - coef[2]


def jeeves_hooke(f, x0, tolerance, coef, r):
    delta = 1.0
    al = 2.0
    x_curr = x0.copy()
    x_best = x0.copy()
    while delta >= tolerance:
        for i in range(2):
            arr = np.array([0, 0], dtype=float)
            arr[i] = delta
            if f(x_curr + arr, coef, r) < f(x_curr, coef, r):
                x_curr += arr
            elif f(x_curr - arr, coef, r) < f(x_curr, coef, r):
                x_curr -= arr
        if f(x_curr, coef, r) < f(x_best, coef, r):
            x_curr, x_best = x_curr + al * (x_curr - x_best), x_curr
        else:
            delta /= al
            x_curr = x_best.copy()
    return x_best


def penalty_method(x, coef, F, f, r=1, c=6, eps=1e-3):
    while F(x, coef, r) - f(x, coef) > eps:
        x = jeeves_hooke(F, x, 1e-9, coef, r)
        r *= c
    return f(x, coef)


def main():
    coef = list(map(float, input().split()))
    x = np.array([2.0, 2.0])
    print(f"{penalty_method(x, coef, F, f):.3f}")


if __name__ == '__main__':
    main()
