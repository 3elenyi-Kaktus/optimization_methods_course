from time import sleep

import numpy as np
# def Nealder_Mead(f, x0, tolerance, coeffs):
#     alpha = 1
#     beta = 0.5
#     gamma = 2
#     while count_tolerance(x0, f, coeffs) >= tolerance:
#         x0 = sorted(x0, key=lambda x: f(x, coeffs))
#         print(x0[0], x0[1], x0[2])
#         x_center = (x0[0] + x0[1]) / 2
#         print(f'Center: {x_center}')
#         x_reflected = x_center + alpha * (x_center - x0[2])
#         print(f'Reflected: {x_reflected}')
#         if f(x_reflected, coeffs) <= f(x0[0], coeffs):
#             x_expanded = x_center + gamma * (x_reflected - x_center)
#             print(f'Expanded: {x_expanded}')
#             if f(x_expanded, coeffs) < f(x0[0], coeffs):
#                 x0[2] = x_expanded
#             else:
#                 x0[2] = x_reflected
#         elif f(x_reflected, coeffs) <= f(x0[1], coeffs):
#             x0[2] = x_reflected
#         else:
#             if f(x_reflected, coeffs) <= f(x0[2], coeffs):
#                 x0[2] = x_reflected
#             x_contraction = x_center + beta * (x_center - x0[2])
#             print(f'Contraction: {x_contraction}')
#             if f(x_contraction, coeffs) <= f(x0[2], coeffs):
#                 x0[2] = x_contraction
#             else:
#                 print('Reduction')
#                 x0[1] = x0[0] + 0.5 * (x0[1] - x0[0])
#                 x0[2] = x0[0] + 0.5 * (x0[2] - x0[0])
#         sleep(0.1)
#         print()
#     x0 = sorted(x0, key=lambda x: f(x, coeffs))
#     return f(x0[0], coeffs)

def f0(x, coef):
    return (4 * (x[0] - coef[0]) ** 2 + (x[1] - coef[1]) ** 2)
def f(x, _):
    return 4*(x[0] - 5)**2 + (x[1] - 6)**2

def f1(x, coef):
    return (x[0] - coef[0]) ** 2 + x[0] * x[1] + coef[1] * (x[1] - 3) ** 2


def count_tolerance(x, f, coeffs):
    x = sorted(x, key=lambda x: f(x, coeffs))
    f_x_center = f((x[0] + x[1]) / 2, coeffs)
    return np.sqrt(((f(x[0], coeffs) - f_x_center) ** 2 + (f(x[1], coeffs) - f_x_center) ** 2 + (f(x[2], coeffs) - f_x_center) ** 2) / 3)


def Nealder_Mead(f, x0, tolerance, coeffs):
    alpha = 1
    beta = 0.5
    gamma = 2
    while count_tolerance(x0, f, coeffs) >= tolerance:
        x0 = sorted(x0, key=lambda x: f(x, coeffs))
        print(x0[0], x0[1], x0[2])
        x_center = (x0[0] + x0[1]) / 2
        print(f'Center: {x_center}')
        x_reflected = x_center + alpha * (x_center - x0[2])
        print(f'Reflected: {x_reflected}')
        if f(x_reflected, coeffs) <= f(x0[0], coeffs):
            x_expanded = x_center + gamma * (x_reflected - x_center)
            print(f'Expanded: {x_expanded}')
            if f(x_expanded, coeffs) < f(x0[0], coeffs):
                x0[2] = x_expanded
            else:
                x0[2] = x_reflected
        elif f(x_reflected, coeffs) <= f(x0[1], coeffs):
            x0[2] = x_reflected
        else:
            if f(x_reflected, coeffs) <= f(x0[2], coeffs):
                x0[2] = x_reflected
            x_contraction = x_center + beta * (x_center - x0[2])
            print(f'Contraction: {x_contraction}')
            if f(x_contraction, coeffs) <= f(x0[2], coeffs):
                x0[2] = x_contraction
            else:
                print('Reduction')
                x0[1] = x0[0] + 0.5 * (x0[1] - x0[0])
                x0[2] = x0[0] + 0.5 * (x0[2] - x0[0])
        sleep(0.1)
        print()
    x0 = sorted(x0, key=lambda x: f(x, coeffs))
    return f(x0[0], coeffs)


def main():
    type = int(input())
    f = f0 if (type == 0) else f1
    coef = [i for i in map(float, input().split())]
    x0 = []
    for k in range(3):
        x0.append(np.array([i for i in map(float, input().split())]))
    x0 = np.array([[8.0, 9.0], [10.0, 11.0], [8.0, 11.0]])
    tol = float(input())
    r1 = Nealder_Mead(f, x0, tol, coef)
    print("{:.10f}".format(r1))


if __name__ == "__main__":
    main()
