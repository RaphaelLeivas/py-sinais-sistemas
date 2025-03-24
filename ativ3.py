import numpy as np
import matplotlib as plt
from sympy import *
from sympy.physics.control.control_plots import pole_zero_plot
from sympy.physics.control.control_plots import bode_plot
from sympy.physics.control.lti import TransferFunction

## Exercicio 1

t, s = symbols('t, s')
funcoes = [exp(3*t) * cos(2*t), t**5 * exp(4*t), sin(5*t) * cos(6*t), sinh(t)]

for f in funcoes:
    print(latex(laplace_transform(f, t, s, noconds=True)))

## Exercicio 2

print("Ex 2 --------------")

funcoes_inv = [1 / (s**2 - 2*s - 3), (5*s + 45) / s**9, -s / (s**3 + 1), (s**3 + 81) / (s**4 + 12 * s**2 + 11)]
for f in funcoes_inv:
    print(latex(inverse_laplace_transform(f, s, t, noconds=True)))

## Exercicio 3

c, C = symbols('c C', cls = Function)

# cond iniciais
y0 = 0
dy_0 = 0

eq_s = Eq(32 * C(s) + 12 * s * C(s) + s**2 * C(s), 2 / s**3)

cs_solucao = solve(eq_s, C(s))[0]

print(latex(cs_solucao))

print(latex(inverse_laplace_transform(cs_solucao, s, t, noconds=True)))

## Exercicio 4

# extrai o numerador e denominador
num, den = fraction(cs_solucao)
ft1 = TransferFunction(num, den, s)
pole_zero_plot(ft1, pole_color="red", grid=False)

## Exercicio 5

bode_plot(ft1, initial_exp=-1, final_exp=3, phase_unit='deg', freq_unit='Hz')