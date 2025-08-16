
# the values of k for which omega(1,k,1) is known

VXXZ_KS = [0, 
    0.321334,

	0.32,
	0.33,
	0.34,
	0.35,
	0.40,
	0.45,
	0.50,

	0.527661,
	0.55,
	0.60,
	0.65,
	0.70,
	0.75,
	0.80,
	0.85,
	0.90,
	0.95,
	1.0,

	1.1,
	1.2,
	1.3,
	1.4,
	1.5,
	1.75,
	2.0,
	2.5,
	3.0,
	4.0,
	5.0
	]

VXXZ_OMEGAS = [2,
	2,

	2.000064,
	2.000100,
	2.000600,
	2.001363,
	2.009541,
	2.023788,
	2.042994,

	2.055322,
	2.066134,
	2.092631,
	2.121734,
	2.153048,
	2.186210,
	2.220929,
	2.256984,
	2.294209,
	2.332440,
	2.371552,
	
	2.452056,
	2.535063,
	2.621644,
	2.708400,
	2.794941,
	3.021591,
	3.250385,
	3.720468,
	4.198809,
	5.171210,
	6.157233
	]

NEW_K_OMEGA_DICT = {
    0.33: 2.000092,
    0.34: 2.000520,
    0.35: 2.001243,
    0.40: 2.009280,
    0.50:  2.042776,
    0.527500: 2.054999,
    0.60: 2.092351,
    0.70: 2.152770,
    0.80: 2.220639,
    0.90: 2.293941,
    1.00: 2.371339,
    1.50: 2.794633,
    2.00: 3.250035
}

def omega(in_a, in_b, in_c, update=True):
    a, b, c = sorted([in_a, in_b, in_c], reverse=True)
    if (c < 0.0001):
        return a + b + c
    if (a - b < 0.0001): 
        return a * rect_omega(c/a, update)
    if (b - c < 0.0001):
        return b * rect_omega(a/b, update)
    
    return min((a-b) + b * rect_omega(c / b, update), (b-c) + c * rect_omega(a / c, update))

def rect_omega(x, update=True):
    ks = VXXZ_KS
    omegas = VXXZ_OMEGAS

    if (len(ks) != len(omegas)):
        raise ValueError('K and Omega Array Lengths Differ')
    
    for i in range(1, len(ks)):
        if ks[i] >= x:
            intLength = ks[i] - ks[i - 1]
            offset = (x - ks[i - 1])/intLength
            prev_omega = NEW_K_OMEGA_DICT.get(ks[i - 1], omegas[i - 1]) if update else omegas[i - 1]
            curr_omega = NEW_K_OMEGA_DICT.get(ks[i], omegas[i]) if update else omegas[i]
            return prev_omega * (1 - offset) + curr_omega * offset
    
    return (x - ks[-1]) + omegas[-1]
