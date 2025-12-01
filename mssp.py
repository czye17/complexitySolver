from functools import partial

import numpy as np

from omega import omega, read_omegas
from scipy.optimize import minimize_scalar

T_RESOLUTION = 1000

def mssp_time(omega_dict, s, t):
    return min(1 + s + t, 2 * (1 - t) + omega(s + t - 1, t, t, omega_dict))

def neg_mssp_time(omega_dict, s, t):
    return - mssp_time(omega_dict, s, t)

def format_res(res, neg=False):
    if neg:
        return -float(res.fun), float(res.x)
    else:
        return float(res.fun), float(res.x)

def mssp_total_time(omega_dict, s):
    res = minimize_scalar(partial(neg_mssp_time, omega_dict, s), bounds=(0, 1))
    return format_res(res, neg=True)
    
def alg_test(omega_dict):
    apsp = omega(1, 1, 1, omega_dict)

    for s in np.linspace(0, 1, endpoint=True, num=100):
        brute_force = min(apsp, 2 + s)
        our_time, _ = mssp_total_time(omega_dict, s)
        opt_time = omega(s, 1, 1, omega_dict)
        if our_time < brute_force:
            print(f'{s:.2f}, {our_time:.4f}, {opt_time:.4f}, {brute_force:.4f}')
    

if __name__ == '__main__':
    omega_dict = read_omegas('/Users/christopherye/Documents/UCSD/Projects/complexitySolver/advxxz.csv')
    alg_test(omega_dict)
