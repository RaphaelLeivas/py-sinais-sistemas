import numpy as np
import matplotlib as plt
from sympy import *
from sympy.physics.control.control_plots import pole_zero_plot
from sympy.physics.control.control_plots import bode_plot
from sympy.physics.control.lti import TransferFunction
from sympy.physics import control

## Exercicio 1
t, s = symbols('t, s')
c, C = symbols('c C', cls = Function)

num1A = s**3 + 2*s**2 + 5*s + 1
den1A = (s + 1) * (s + 2) * (s + 3) * (s+4)
ft1A = TransferFunction(num1A, den1A, s)

ht1A = inverse_laplace_transform(num1A / den1A, s, t, noconds=True)
print(latex(ht1A))

BUG()
