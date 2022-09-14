from sympy import symbols, I, solve, Matrix
from pickle import load as pickle_load
from src.modules.member import member as member_type

c1, d1 = symbols("c_traveling, d_traveling")
c2, d2 = symbols("c_evanescent, d_evanescent")
rho1, Area1, E1, I1 = symbols("rho1, Area1, E1, I1")


def _subs(eqns: list, m1: member_type) -> list:
    subs_dict = {
        rho1: m1.rho,
        Area1: m1.cross_section_area,
        E1: m1.youngs_modulus,
        I1: m1.inertia,
    }
    for i in range(len(eqns)):
        for j in range(len(eqns[i])):
            eqns[i][j] = eqns[i][j].subs(subs_dict)
    return eqns


def _get_soln(eqns):
    travelling_eqs, evanscent_eqs = eqns

    soln_1 = solve(travelling_eqs, (c1, d1))
    c1_sol = soln_1[c1]
    d1_sol = soln_1[d1]

    soln_2 = solve(evanscent_eqs, (c2, d2))
    c2_sol = soln_2[c2]
    d2_sol = soln_2[d2]

    reflection_matrix = Matrix([[c1_sol, d1_sol], [c2_sol, d2_sol]])
    return reflection_matrix


def get_r_of_free_end(m1: member_type) -> tuple:
    file = open("src/equations/free_end", mode="rb")
    eqns = pickle_load(file)
    eqns = _subs(eqns, m1)
    return _get_soln(eqns)


def get_r_of_fixed_end(m1: member_type) -> tuple:
    file = open("src/equations/fixed_end", mode="rb")
    eqns = pickle_load(file)
    eqns = _subs(eqns, m1)
    return _get_soln(eqns)
