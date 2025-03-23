from sympy import *

## Exercicio 1

t = symbols('t')
y = Function('y')(t)
x = Function('x')(t)

QN_coeffs = [9, 6, 1]
PN_coeffs = [7, 2]
cond_iniciais = [3, -1]

# define a EDO
edo_y = sum(coeff * Derivative(y, t, n) for n, coeff in enumerate(QN_coeffs))
edo_x = sum(coeff * Derivative(x, t, n) for n, coeff in enumerate(PN_coeffs))
edo_completa = Eq(edo_y, edo_x)
print(latex(edo_completa))

eq_homog = Eq(edo_y, 0) #x(t) = 0, entrada zero

# edo linear homogenea de segunda ordem: resolve com polinomio característico
lambda_ = symbols('lambda')
eq_carac = sum(coeff * lambda_**n for n, coeff in enumerate(QN_coeffs))
print(latex(eq_carac))
raizes = solve(eq_carac, lambda_)
print(raizes)

# analisa as raizes da eq caracteristica para montar a solucao geral
sol_geral = 0
constantes = []
c_counter = 1

raizes_mult = roots(eq_carac, lambda_)
for raiz, mult in raizes_mult.items():
    for m in range(mult):
        constante = symbols(f'c{c_counter}')
        constantes.append(constante)
        sol_geral += constante * t**m * exp(raiz * t)
        c_counter += 1


# cond iniciais
condicoes = []
for ordem, valor in enumerate(cond_iniciais):
    condicoes.append(Eq(sol_geral.diff(t, ordem).subs(t, 0), valor))

sistema_equacoes = []
for cond in condicoes:
    sistema_equacoes.append(cond.lhs - cond.rhs)

solucao_sistema = solve(sistema_equacoes, constantes)

solucao_final = sol_geral.subs(solucao_sistema)

print(latex(solucao_final))

## Exercicio 2

def resposta_ao_impulso(Ps, Qs):
    t = symbols('t', real=True, positive=True)
    s = symbols('s')

    Hs = Ps / Qs

    h_t = inverse_laplace_transform(Hs, s, t)

    return h_t

s = symbols('s')

QN_coeffs = [-4, 3, 1]
PN_coeffs = [-2, 0, 1]
edo_y = sum(coeff * Derivative(y, t, n) for n, coeff in enumerate(QN_coeffs))
edo_x = sum(coeff * Derivative(x, t, n) for n, coeff in enumerate(PN_coeffs))
edo_completa = Eq(edo_y, edo_x)
print(latex(edo_completa))

Qs = -4 + 3*s + s**2
Ps = -2 + s**2

h_t = resposta_ao_impulso(Ps, Qs)

# verifica se é instantaneo
grau_p = degree(Ps, s)
grau_q = degree(Qs, s)

if grau_p == grau_q:
    b_0 = LC(Ps, s)
    a_0 = LC(Qs, s)

    termo_dirac = b_0 / a_0
    print(latex(termo_dirac + h_t))
else:
    print(latex(h_t))

## Exercicio 3

t, tau = symbols('t, tau')

h_tau = 3*exp(-6 * tau) + exp(-tau)
x_tau = exp(tau)

convolution = integrate(x_tau * h_tau.subs(tau, t - tau), (tau, 0, t))
print(latex(convolution))

t_value = 1
y_value = convolution.subs(t, t_value)

print(f"y({t_value}) = ", y_value.evalf())