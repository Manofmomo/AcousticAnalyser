from sympy import symbols, I, solve, Matrix
from pickle import load as pickle_load
from src.modules.member import member as member_type
from typing import List, Tuple

a1, b1, c1, d1 = symbols("a_traveling,  b_traveling, c_traveling, d_traveling")
a2, b2, c2, d2 = symbols("a_evanescent, b_evanescent, c_evanescent, d_evanescent")
rho1, Area1, E1, I1 = symbols("rho1, Area1, E1, I1")
rho2, Area2, E2, I2 = symbols("rho2, Area2, E2, I2")
theta = symbols("theta")


def _subs(eqns: List[Matrix], m1: member_type, m2: member_type, theta: float) -> list:
    subs_dict = {
        rho1: m1.rho,
        Area1: m1.cross_section_area,
        E1: m1.youngs_modulus,
        I1: m1.inertia,
        rho2: m2.rho,
        Area2: m2.cross_section_area,
        E2: m2.youngs_modulus,
        I2: m2.inertia,
        theta: theta,
    }
    for i in range(len(eqns)):
        eqns[i] = eqns[i].subs(subs_dict)
    return eqns


def _get_soln(eqns: List[Matrix]) -> Tuple(Matrix):
    M1, M2, M3, M4, M5, M6 = eqns

    transmission_matrix = (M5.inv() * M4 - M2.inv() * M1).inv() * (
        M5.inv() * M6 - M2.inv() * M3
    )
    reflection_matrix = (M1.inv() * M2 - M4.inv() * M5).inv() * (
        M4.inv() * M6 - M1.inv() * M3
    )

    return (reflection_matrix, transmission_matrix)


def get_rt_of_cross_section(m1: member_type, m2: member_type, theta: float) -> tuple:
    file = open("src/equations/cross_section", mode="rb")
    eqns = pickle_load(file)
    eqns = _subs(eqns, m1, m2, theta)
    return _get_soln(eqns)
