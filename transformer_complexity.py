from omega import omega, read_omegas

def main():
    omega_dict = read_omegas('/Users/christopherye/Documents/UCSD/Projects/complexitySolver/advxxz.csv')

    sum_mm_time = omega(1, 2, 1, omega_dict)
    print(sum_mm_time)
    each_mm_time = 1 + omega(1, 1, 1, omega_dict)
    print(each_mm_time)

    T = 100
    for x in range(int(0.3*T), T):
        a = x/T
        rect_mm = omega(1, x/T, 1, omega_dict)
        worse_rect_mm1 = max(2,omega(1, a, 1, omega_dict) - (1-a))
        worse_rect_mm2 = max(2,omega(1, a, a, omega_dict))
        worse_rect_mm = max(worse_rect_mm1, worse_rect_mm2)
        which_worst = worse_rect_mm2 > worse_rect_mm1
        print(f'{a:3f} {rect_mm:.6f} {worse_rect_mm:.6f} {rect_mm - worse_rect_mm:6f} {which_worst}')


if __name__ == '__main__':
    main()