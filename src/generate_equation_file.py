#%%
import sympy
from sympy.parsing.mathematica import parse_mathematica as parse
from sympy import I, sin, cos, tan
import pickle
from typing import Tuple
import logging
#%%
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s %(name)s %(levelname)s:%(message)s"
)
logger = logging.getLogger(__name__)
# This initial code is developed using Euler-Bernoulli beam mode


def get_equation_free_end() -> sympy.Matrix:
    reflection = sympy.Matrix([[-I, -1 + I, 0], [-1 + I, I, 0], [0, 0, -1]])

    return reflection


def get_equation_fixed_end() -> sympy.Matrix:
    reflection = sympy.Matrix([[-I, -1 + I, 0], [1 - I, I, 0], [0, 0, 1]])

    return reflection


def get_equation_cross() -> Tuple[sympy.Matrix]:
    """Gets the M0-M6 of a joint assuming Wave equation to be A*e^(-I*k*Zeta) + B*e^(-k*Zeta) + C*e^(I*k*Zeta) + D*e^(k*Zeta)"""
    #%%
    w = sympy.symbols("w")
    density1, area1, E1, I1, H1 = sympy.symbols("density1, area1, E1, I1, H1")
    theta = sympy.symbols("theta")

    # tO REMOVE
    pi=3.1416
    theta=90*pi/180
    E1=206e9
    density1=7800
    H1=1.27e-2        
    I1=(H1**4)/12
    area1=H1*H1
    w=2*pi*20
    # Till here

    C = (E1 / density1) ** 0.5
    K1 = (I1 / area1) ** 0.5

    Zeta = H1 / K1
    omega = (w * K1 / C) ** 0.5
    
    #%%
    M1 = sympy.Matrix(parse("{{I*Sin[theta], -Sin[theta], ((-I)*Cos[theta])/omega}, {(-I)*Cos[theta], Cos[theta], ((-I)*Sin[theta])/omega}, {(-omega^2 - (I/2)*Zeta*omega^3*Tan[theta/2])/omega^2, (omega^2 + (Zeta*omega^3*Tan[theta/2])/2)/omega^2, 0}}"))
    M2 = sympy.Matrix(parse("{{0, 0, -(((-I)*omega^2 + Zeta*omega^4)/omega^3)}, {-(((-I)*omega^3 + Zeta*omega^4 + (I/2)*Zeta^2*omega^5*Tan[theta/2])/omega^3), -((omega^3 + Zeta*omega^4 + (Zeta^2*omega^5*Tan[theta/2])/2)/omega^3), 0}, {-((omega^2 + (I/6)*Zeta^3*omega^5 + (I/2)*Zeta*omega^3*Tan[theta/2])/omega^2), -((-omega^2 + (Zeta^3*omega^5)/6 - (Zeta*omega^3*Tan[theta/2])/2)/omega^2), 0}}"))
    M3 = sympy.Matrix(parse("{{0, 0, -((I*omega^2 + Zeta*omega^4)/omega^3)}, {-((I*omega^3 + Zeta*omega^4 - (I/2)*Zeta^2*omega^5*Tan[theta/2])/omega^3), -((-omega^3 + Zeta*omega^4 - (Zeta^2*omega^5*Tan[theta/2])/2)/omega^3), 0}, {-((omega^2 - (I/6)*Zeta^3*omega^5 - (I/2)*Zeta*omega^3*Tan[theta/2])/omega^2), -((-omega^2 - (Zeta^3*omega^5)/6 + (Zeta*omega^3*Tan[theta/2])/2)/omega^2), 0}}"))
    M4 = sympy.Matrix(parse("{{0, 0, 1}, {1, 1, 0}, {-I, -1, 0}}"))
    M5 = sympy.Matrix(parse("{{Sin[theta] + (I/2)*Zeta*omega*Sin[theta]*Tan[theta/2], Sin[theta] + (Zeta*omega*Sin[theta]*Tan[theta/2])/2, Cos[theta]}, {Cos[theta] + (I/2)*Zeta*omega*Tan[theta/2] + (I/2)*Zeta*omega*Cos[theta]*Tan[theta/2], Cos[theta] + (Zeta*omega*Tan[theta/2])/2 + (Zeta*omega*Cos[theta]*Tan[theta/2])/2, -Sin[theta]}, {I, 1, 0}}"))
    M6 = sympy.Matrix(parse("{{Sin[theta] - (I/2)*Zeta*omega*Sin[theta]*Tan[theta/2], Sin[theta] - (Zeta*omega*Sin[theta]*Tan[theta/2])/2, Cos[theta]}, {Cos[theta] - (I/2)*Zeta*omega*Tan[theta/2] - (I/2)*Zeta*omega*Cos[theta]*Tan[theta/2], Cos[theta] - (Zeta*omega*Tan[theta/2])/2 - (Zeta*omega*Cos[theta]*Tan[theta/2])/2, -Sin[theta]}, {-I, -1, 0}}"))

    N1 = sympy.Matrix(parse("{{I*Sin[theta], -Sin[theta], ((-I)*Cos[theta])/omega}, {(-I)*Cos[theta], Cos[theta], ((-I)*Sin[theta])/omega}, {(-omega^2 - (I/2)*Zeta*omega^3*Tan[theta/2])/omega^2, (omega^2 + (Zeta*omega^3*Tan[theta/2])/2)/omega^2, 0}}"))
    N2 = sympy.Matrix(parse("{{0, 0, -(((-I)*omega^2 + Zeta*omega^4)/omega^3)}, {-(((-I)*omega^3 + Zeta*omega^4 + (I/2)*Zeta^2*omega^5*Tan[theta/2])/omega^3), -((omega^3 + Zeta*omega^4 + (Zeta^2*omega^5*Tan[theta/2])/2)/omega^3), 0}, {-((omega^2 + (I/6)*Zeta^3*omega^5 + (I/2)*Zeta*omega^3*Tan[theta/2])/omega^2), -((-omega^2 + (Zeta^3*omega^5)/6 - (Zeta*omega^3*Tan[theta/2])/2)/omega^2), 0}}"))
    N3 = sympy.Matrix(parse("{{0, 0, -((I*omega^2 + Zeta*omega^4)/omega^3)}, {-((I*omega^3 + Zeta*omega^4 - (I/2)*Zeta^2*omega^5*Tan[theta/2])/omega^3), -((-omega^3 + Zeta*omega^4 - (Zeta^2*omega^5*Tan[theta/2])/2)/omega^3), 0}, {-((omega^2 - (I/6)*Zeta^3*omega^5 - (I/2)*Zeta*omega^3*Tan[theta/2])/omega^2), -((-omega^2 - (Zeta^3*omega^5)/6 + (Zeta*omega^3*Tan[theta/2])/2)/omega^2), 0}}"))
    N4 = sympy.Matrix(parse("{{0, 0, 1}, {1, 1, 0}, {-I, -1, 0}}"))
    N5 = sympy.Matrix(parse("{{Sin[theta] + (I/2)*Zeta*omega*Sin[theta]*Tan[theta/2], Sin[theta] + (Zeta*omega*Sin[theta]*Tan[theta/2])/2, Cos[theta]}, {Cos[theta] + (I/2)*Zeta*omega*Tan[theta/2] + (I/2)*Zeta*omega*Cos[theta]*Tan[theta/2], Cos[theta] + (Zeta*omega*Tan[theta/2])/2 + (Zeta*omega*Cos[theta]*Tan[theta/2])/2, -Sin[theta]}, {I, 1, 0}}"))
    N6 = sympy.Matrix(parse("{{Sin[theta] - (I/2)*Zeta*omega*Sin[theta]*Tan[theta/2], Sin[theta] - (Zeta*omega*Sin[theta]*Tan[theta/2])/2, Cos[theta]}, {Cos[theta] - (I/2)*Zeta*omega*Tan[theta/2] - (I/2)*Zeta*omega*Cos[theta]*Tan[theta/2], Cos[theta] - (Zeta*omega*Tan[theta/2])/2 - (Zeta*omega*Cos[theta]*Tan[theta/2])/2, -Sin[theta]}, {-I, -1, 0}}"))
    #%%
    return [M1, M2, M3, M4, M5, M6, N1, N2, N3, N4, N5, N6]


if __name__ == "__main__":
    obj = get_equation_cross()
    filehandler = open("equations/cross_section", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    logger.debug("cross_section generated")

    obj = get_equation_free_end()
    filehandler = open("equations/free_end", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    logger.debug("free_end generated")

    obj = get_equation_fixed_end()
    filehandler = open("equations/fixed_end", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    logger.debug("fixed_end generated")

# %%
