import numpy as np
import scipy.optimize as scopt

def main():
    c = -np.array([4., 1., -4., -4.])
    A_ub = np.array([[3., 0., 5., 2.],
                     [1., 2., -1., 0.],])
    b_ub = np.array([3., 1.])
    # проставляем индикаторы для типа переменных, 0 - непрерывное, 1 - целое число
    var_types = [1, 1, 1, 1]
    # также указываем границы, в том числе и для целочисленных переменных
    bounds = [(0, None), (0, None), (0, None), (0, None)]
    res_milp = scopt.linprog(c=c, A_eq=A_ub, b_eq=b_ub, bounds=bounds, integrality=var_types, method='highs')
    print(f"Решение: x = {list(np.round(res_milp['x'], 2))}, f = {-res_milp['fun']}, {res_milp['message']}")


if __name__ == '__main__':
    main()
