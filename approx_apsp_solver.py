from functools import partial

import numpy as np
from scipy.optimize import minimize_scalar

from omega import omega, read_omegas

X_RESOLUTION = 1000000
L_RESOLUTION = 200
USE_NEW_VALUES = True

def opt_x_numeric(func):
    opt_x = -1
    opt_time = 3
    for x in np.linspace(0, 1, X_RESOLUTION):
        time = func(x)
        if time < opt_time:
            opt_time = time
            opt_x = x
    return float(opt_time), float(opt_x)

def format_res(res):
    return float(res.fun), float(res.x)

def mult_approx_opt(omega_dict):
    # return minimize(mult_approx_time, 0.5, method = 'Nelder-Mead')
    res = minimize_scalar(partial(mult_approx_time, omega_dict))
    return format_res(res)
    
def mult_approx_time(omega_dict, x):
    mm_time = omega(1, 1-x, 1, omega_dict)
    time = max(mm_time, 1.5 + x)
    return time

def det_mult_approx_opt(omega_dict, k=7):
    res = minimize_scalar(partial(det_mult_approx_time, k, omega_dict))
    return format_res(res)

def det_mult_approx_time(k, omega_dict, x):
    mm_time = omega(1, 1-x, 1, omega_dict)
    time = max(mm_time, 1 + 2*x, 2 + x/k, 1 + ((3*k - 4)*x)/k)
    return time

def add_approx_opt(omega_dict, k=4, bounded=False):
    res = minimize_scalar(partial(add_approx_time, k, omega_dict, bounded))
    return format_res(res)

def add_approx_time(k, omega_dict, bounded, x):
    in_a = 1-((k-2)*x)/(k+2)
    in_b = 1 - x
    in_c = 1 - ((k-4)*x)/(k+2)
    if bounded:
        mm_time = omega(in_a, in_b, in_c, omega_dict)
    else:
        in_mu = 1
        mm_time = ((in_a + in_b + in_mu) + omega(in_a, in_b, in_c, omega_dict))/2

    time = max(mm_time, 2 + (2*x)/(k+2))
    return time

def weighted_add_approx_opt(omega_dict, k=2):
    res = minimize_scalar(partial(weighted_add_approx_time, k, omega_dict))
    return format_res(res)

def weighted_add_approx_time(k, omega_dict, x):
    if k == 1:
        mm_time = omega(1, 1-x, 1, omega_dict)
        time = max(mm_time, 2 + x/2)
        return time 
    else:
        in_a = 1-((k-1)*x)/(k+1)
        in_b = 1 - x
        in_c = 1 - ((k-2)*x)/(k+1)
        mm_time = omega(in_a, in_b, in_c, omega_dict)

        time = max(mm_time, 2 + x/(k+1))
        return time

def mult_approx_long_opt(omega_dict, k=4):
    opt_x, opt_l = -1, -1
    opt_time = 3
    for l in range(20, L_RESOLUTION):
        res = minimize_scalar(partial(mult_approx_long_time, omega_dict, k, l))
        if res.fun < opt_time:
            opt_time = res.fun
            opt_x = res.x
            opt_l = l
    return float(opt_time), float(opt_x), opt_l

def mult_approx_long_time(omega_dict, k, l, x):
    in_a = 1 - ((k - 2)*x)/(2*(l - k + 2))
    in_b = 1 - x
    in_c = 1 - ((k - 4)*x)/(2*(l - k + 2))
    mm_time = omega(in_a, in_b, in_c, omega_dict)
    time = max(mm_time, 2 + x/(l - k + 2), 1.5 + x)
    return time

def generate_results():
    omega_dict = read_omegas('/Users/christopherye/Documents/UCSD/Projects/complexitySolver/advxxz.csv')

    print('\n\n---- 2 Multiplicative Approximation ----')
    print(mult_approx_opt(omega_dict))

    print('\n\n---- Deterministic 2 Multiplicative Approximation ----')
    print(det_mult_approx_opt(omega_dict))

    print('\n\n---- Additive Approximation----')
    for k in range(4, 12, 2):
        print(add_approx_opt(omega_dict, k))

    print('\n\n---- Bounded Additive Approximation----')
    for k in range(4, 14, 2):
        print(add_approx_opt(omega_dict, k, bounded=True))

    print('\n\n---- Weighted Additive Approximation----')
    for k in range(1, 6):
        print(weighted_add_approx_opt(omega_dict, k))

    print('\n\n---- Long Path Multiplicative Approximation----')
    for k in range(4, 14, 2):
        print(mult_approx_long_opt(omega_dict, k))
    print('\n\n')


if __name__ == '__main__':
    generate_results()
