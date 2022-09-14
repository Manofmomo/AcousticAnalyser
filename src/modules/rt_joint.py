from sympy import symbols, I,solve, Matrix
from pickle import load as pickle_load
from src.modules.member import member as member_type

a1, b1, c1, d1= symbols("a_traveling,  b_traveling, c_traveling, d_traveling")
a2, b2, c2, d2= symbols("a_evanescent, b_evanescent, c_evanescent, d_evanescent")
rho1, Area1, E1, I1 = symbols("rho1, Area1, E1, I1")
rho2, Area2, E2, I2 = symbols("rho2, Area2, E2, I2")

def _subs(eqns: list, m1: member_type, m2: member_type) -> list:
    subs_dict = {rho1:m1.rho,Area1:m1.cross_section_area,E1:m1.youngs_modulus,I1:m1.inertia,rho2:m2.rho, Area2:m2.cross_section_area, E2:m2.youngs_modulus, I2:m2.inertia}
    for i in range(len(eqns)):
        for j in range(len(eqns[i])):
            eqns[i][j] = eqns[i][j].subs(subs_dict)
    return eqns

def _get_soln(eqns):
    travelling_eqs, evanscent_eqs = eqns

    soln_1 = solve(travelling_eqs,(a1, b1, c1, d1))
    a1_sol = soln_1[a1]
    b1_sol = soln_1[b1]
    c1_sol = soln_1[c1]
    d1_sol = soln_1[d1]

    soln_2 = solve(evanscent_eqs,(a2, b2, c2, d2))
    a2_sol = soln_2[a2]
    b2_sol = soln_2[b2]
    c2_sol = soln_2[c2]
    d2_sol = soln_2[d2]

    transmission_matrix = Matrix([[a1_sol,b1_sol],[a2_sol, b2_sol]])
    reflection_matrix = Matrix([[c1_sol,d1_sol],[c2_sol, d2_sol]])
    return (reflection_matrix, transmission_matrix)

def get_rt_of_cross_section(m1: member_type, m2: member_type, theta: float) -> tuple:
    file = open('src/equations/cross_section',mode='rb')
    eqns = pickle_load(file)
    eqns = _subs(eqns,m1,m2)
    return _get_soln(eqns)