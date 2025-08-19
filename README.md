# complexitySolver

Script to solve complexity optimization problems.

For complexities involving matrix multiplication, `omega(a, b, c)` denotes the running time to compute the product of a $n^{a} \times n^{b}$ matrix and a $n^{b} \times n^{c}$ matrix.
Assuming $a > b > c$, the function either computes $n^{a - b}$ rectangular products of shape $n^{\omega(b, b, c)}$ or $n^{b - c}$ rectangular products of shape $n^{\omega(a, c, c)}$.

Code is adapted from https://www.ocf.berkeley.edu/~vdbrand/complexity/ by Jan van den Brand and updated to include most recent bounds on fast matrix multiplication (`More Asymmetry Yields Faster Matrix Multiplication` by Alman, Duan, Vassilevska-Williams, Xu, Xu, Zhou SODA 2025)

Packages required: numpy and scipy.

Run `python approx_apsp_solver.py`
