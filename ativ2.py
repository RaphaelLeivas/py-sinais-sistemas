from sympy import *

t = symbols('t')
y = Function('y')(t)

QN_coeffs = [9, 6, 1]
cond_iniciais = [5, -2]

# define a edo
edo = sum(coeff * Derivative(y, t, n) for n, coeff in enumerate(QN_coeffs))
eq_homog = Eq(edo, 0) #x(t) = 0, entrada zero


# edo linear homogenea de segunda ordem: resolve com polinomio característico
lambda_ = symbols('lambda')
eq_carac = sum(coeff * lambda_**n for n, coeff in enumerate(QN_coeffs))
raizes = solve(eq_carac, lambda_)

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

