import numpy as np

from omega import omega, read_omegas

def main():
    omega_dict = read_omegas('/Users/christopherye/Documents/UCSD/Projects/complexitySolver/advxxz.csv')
    apsp = omega(1, 1, 1, omega_dict)

    for s in np.linspace(0, 1):
        brute_force = min(apsp, 2 + s)
        our_time = 1 + 0.5 * (s + omega(1, 1, s, omega_dict))
        if our_time < brute_force:
            print(f's: {s:.2f} gap: {brute_force - our_time:.6f}')

if __name__ == '__main__':
    main()
