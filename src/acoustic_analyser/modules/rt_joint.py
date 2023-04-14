from sympy import symbols, I, Matrix
from sympy.core.parameters import evaluate
from pickle import load as pickle_load
from acoustic_analyser.modules.member import member as member_type
from typing import List, Tuple
import logging
import os

logger = logging.getLogger("acoustic_analyser")


a1, b1, c1, d1 = symbols("a_traveling,  b_traveling, c_traveling, d_traveling")
a2, b2, c2, d2 = symbols("a_evanescent, b_evanescent, c_evanescent, d_evanescent")
density1, area1, E1, I1, L1, H1 = symbols("density1, area1, E1, I1, L1, H1")
density2, area2, E2, I2, L2, H2 = symbols("density2, area2, E2, I2, L2, H2")
theta_sym = symbols("theta")
w_sym = symbols("w")

base_path = os.path.join(os.path.dirname(__file__), "..", "equations")
file = open(os.path.join(base_path, "two_member.pkl"), mode="rb")
eqns = pickle_load(file)
logger.debug("Joint Equation Files Loaded")
file.close()


def _subs(
    eqns: List[Matrix], m1: member_type, m2: member_type, theta: float, w: float
) -> list:
    subs_dict = {
        density1: m1.density,
        area1: m1.cross_section_area,
        E1: m1.youngs_modulus,
        I1: m1.inertia,
        L1: m1.length,
        H1: m1.height,
        density2: m2.density,
        area2: m2.cross_section_area,
        E2: m2.youngs_modulus,
        I2: m2.inertia,
        L2: m2.length,
        H2: m2.height,
        theta_sym: theta,
        w_sym: w,
    }
    eqns_subs: List[Matrix] = [None] * len(eqns)
    for i in range(len(eqns)):
        eqns_subs[i] = eqns[i].subs(subs_dict)
        eqns_subs[i].simplify()
    return eqns_subs


def _get_soln(eqns: List[Matrix]) -> Tuple[Matrix]:
    M1, M2, M3, M4, M5, M6, N1, N2, N3, N4, N5, N6 = eqns
    logger.debug("M1 :")
    logger.debug(M1)
    transmission_matrix_12 = (N5.inv() * N4 - N2.inv() * N1).inv() * (
        N5.inv() * N6 - N2.inv() * N3
    )
    logger.debug("transmission_matrix_12 Calculated")
    reflection_matrix_11 = (N1.inv() * N2 - N4.inv() * N5).inv() * (
        N4.inv() * N6 - N1.inv() * N3
    )
    logger.debug("reflection_matrix_11 Calculated")
    transmission_matrix_21 = (M5.inv() * M4 - M2.inv() * M1).inv() * (
        M5.inv() * M6 - M2.inv() * M3
    )
    logger.debug("transmission_matrix_21 Calculated")
    reflection_matrix_22 = (M1.inv() * M2 - M4.inv() * M5).inv() * (
        M4.inv() * M6 - M1.inv() * M3
    )
    logger.debug("reflection_matrix_22 Calculated")

    return (
        reflection_matrix_11.evalf(),
        reflection_matrix_22.evalf(),
        transmission_matrix_12.evalf(),
        transmission_matrix_21.evalf(),
    )


def get_rt_of_two_member(
    m1: member_type, m2: member_type, theta: float, w: float
) -> tuple:
    global eqns
    eqns_subs = _subs(eqns, m1, m2, theta, w)
    logger.debug("M0-M6 N0-N6 Substituted")
    reflection_transmission = _get_soln(eqns_subs)
    logger.debug("Reflection Transmission Calculated")
    return reflection_transmission
