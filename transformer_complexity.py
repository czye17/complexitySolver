from omega import omega, read_omegas


def main():
    omega_dict = read_omegas('/Users/christopherye/Documents/UCSD/Projects/complexitySolver/advxxz.csv')

    sum_mm_time = omega(1, 2, 1, omega_dict)
    print(sum_mm_time)
    each_mm_time = 1 + omega(1, 1, 1, omega_dict)
    print(each_mm_time)


if __name__ == '__main__':
    main()