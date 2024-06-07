def f0(x, coef):
    return coef[0] * x ** 2 + coef[1] * x + coef[2]


def f1(x, coef):
    return coef[0] * x ** 4 + coef[1] * x ** 3 + coef[2] * x ** 2 + coef[3] * x + coef[4]


def fib(n):
    if n in [0, 1]:
        return 1
    return fib(n - 1) + fib(n - 2)


def fib_search(f, bounds, tol, coef, max_eps=0.01):
    fib_barrier = (bounds[1] - bounds[0]) / tol
    iter = 0
    # print(fib_barrier)
    while fib(iter) < fib_barrier:
        # print(iter, fib(iter))
        iter += 1
    # print(iter, fib(iter))

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
        # print(f'Iter: {k}\nLeft: {left}\nRight: {right}')
        k += 1

    right = left + max_eps
    f_right = f(right, coef)
    if f_left >= f_right:
        bounds[0] = left
    else:
        bounds[1] = right

    # print(f'\nPoint: {(bounds[0] + bounds[1]) / 2}\nF: {f((bounds[0] + bounds[1]) / 2, coef)}')
    return (bounds[0] + bounds[1]) / 2


def tests():
    type = 1
    f = f0 if (type == 0) else f1
    coef = [i for i in map(float, [2, 10, 17, -25, -24])]
    bounds = [0, 0]
    bounds[0], bounds[1], tol = map(float, [-2, -1, 0.0007])
    r1 = fib_search(f, bounds, tol, coef)
    print("{:.10f}".format(r1))

    type = 0
    f = f0 if (type == 0) else f1
    coef = [i for i in map(float, [23, 39, -34])]
    bounds = [0, 0]
    bounds[0], bounds[1], tol = map(float, [-2, -1, 0.0004])
    r1 = fib_search(f, bounds, tol, coef)
    print("{:.10f}".format(r1))

    type = 0
    f = f0 if (type == 0) else f1
    coef = [i for i in map(float, [41, 45, -32])]
    bounds = [0, 0]
    bounds[0], bounds[1], tol = map(float, [-1, 0, 0.0002])
    r1 = fib_search(f, bounds, tol, coef)
    print("{:.10f}".format(r1))


if __name__ == "__main__":
    tests()
    type = int(input())
    f = f0 if (type == 0) else f1
    coef = [i for i in map(float, input().split())]
    bounds = [0, 0]
    bounds[0], bounds[1], tol = map(float, input().split())
    r1 = fib_search(f, bounds, tol, coef)
    print("{:.10f}".format(r1))
