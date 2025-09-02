import csv

def read_omegas(path):
    with open(path, mode='r') as infile:
        reader = csv.reader(infile)
        omega_dict = {float(rows[0]): float(rows[1]) for rows in reader}
        return omega_dict

def omega(in_a, in_b, in_c, omega_dict):
    a, b, c = sorted([in_a, in_b, in_c], reverse=True)
    if (c < 0.0001):
        return a + b + c
    if (a - b < 0.0001): 
        return a * rect_omega(c/a, omega_dict)
    if (b - c < 0.0001):
        return b * rect_omega(a/b, omega_dict)
    
    return min((a-b) + b * rect_omega(c / b, omega_dict), (b-c) + c * rect_omega(a / c, omega_dict))

def rect_omega(x, omega_dict):
    i = 0
    for k, om in omega_dict.items():
        if i >= 1 and k >= x:
            intLength = k - prev_k
            offset = (x - prev_k)/intLength
            return prev_om * (1 - offset) + om * offset
        
        prev_k = k
        prev_om = om
        i += 1
    
    return (x - k) + om
