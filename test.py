from sympy import *

x, t, z, nu = symbols('x t z nu')
init_printing(use_unicode=True)
res = diff(sin(x)*exp(x), x)

res = integrate(sin(x**2), (x, -oo, oo))

print(solve(x**2 - 2, x))