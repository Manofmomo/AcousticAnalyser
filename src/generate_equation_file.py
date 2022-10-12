from locale import MON_1
from zipfile import ZIP_DEFLATED
import sympy
from sympy import I
import pickle

# This initial code is developed using Euler-Bernoulli beam mode


def get_equation_free_end():

    x = sympy.symbols("x")
    w = sympy.symbols("w")
    density1, area1, E1, I1 = sympy.symbols("density1, area1, E1, I1")

    K1 = (density1 * area1 / (E1 * I1)) ** (0.25) * (w**0.5)

    c1, d1 = sympy.symbols("c_traveling, d_traveling")
    c2, d2 = sympy.symbols("c_evanescent, d_evanescent")
    params = [(c1, d1), (c2, d2)]

    travelling_wave = sympy.exp(-I * K1 * x)
    evanescent_wave = sympy.exp(-1 * K1 * x)
    waves = [travelling_wave, evanescent_wave]

    final_equations = []

    for i in range(2):
        incoming_wave = waves[i]
        c, d = params[i]

        v = incoming_wave + c * sympy.exp(I * K1 * x) + d * sympy.exp(K1 * x)
        vdiffs = [None] * 4
        vdiffs[0] = v
        for i in range(3):
            vdiffs[i + 1] = sympy.diff(vdiffs[i], x)

        eq1 = vdiffs[2].subs(x, 0)
        eq2 = vdiffs[3].subs(x, 0)

        final_equations.append([eq1, eq2])

    return final_equations


def get_equation_fixed_end():

    x = sympy.symbols("x")
    w = sympy.symbols("w")
    density1, area1, E1, I1 = sympy.symbols("density1, area1, E1, I1")

    K1 = (density1 * area1 / (E1 * I1)) ** (0.25) * (w**0.5)

    c1, d1 = sympy.symbols("c_traveling, d_traveling")
    c2, d2 = sympy.symbols("c_evanescent, d_evanescent")
    params = [(c1, d1), (c2, d2)]

    travelling_wave = sympy.exp(-I * K1 * x)
    evanescent_wave = sympy.exp(-1 * K1 * x)
    waves = [travelling_wave, evanescent_wave]

    final_equations = []

    for i in range(2):
        incoming_wave = waves[i]
        c, d = params[i]

        v = incoming_wave + c * sympy.exp(I * K1 * x) + d * sympy.exp(K1 * x)
        vdiffs = [None] * 4
        vdiffs[0] = v
        for i in range(3):
            vdiffs[i + 1] = sympy.diff(vdiffs[i], x)

        eq1 = vdiffs[0].subs(x, 0)
        eq2 = vdiffs[1].subs(x, 0)

        final_equations.append([eq1, eq2])

    return final_equations


# Assuming Wave equation to be A*e^(-I*k*x) + B*e^(-k*x) + C*e^(I*k*x) + D*e^(k*x)
def get_equation_cross():

    w = sympy.symbols("w")
    density1, area1, E1, I1, L1 = sympy.symbols("density1, area1, E1, I1, L1")

    H1 = 12 * I1 / (L1**3)
    C = (E1 / density1) ** 0.5
    K1 = (L1 / area1) ** 0.5

    zeta = H1 / K1
    omega = (w * K1 / C) ** 0.5
    theta = sympy.symbols("theta")

    M1 = sympy.Matrix(
        [
            [0, 0, zeta * omega - I / omega],
            [zeta * omega - I, zeta * omega + 1, 0],
            [
                1 + 0.5 * I * zeta * sympy.tan(theta / 2),
                -1 - 0.5 * I * zeta * sympy.tan(theta / 2),
                0,
            ],
        ]
    )

    M2 = sympy.Matrix(
        [
            [-I * sympy.sin(theta), sympy.sin(theta), I * sympy.cos(theta) / omega],
            [
                0.5
                * I
                * (
                    (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                    + 2 * sympy.cos(theta)
                ),
                0.5 * (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                - sympy.cos(theta),
                I * sympy.sin(theta) / omega,
            ],
            [
                (1 / 6)
                * (
                    I * (zeta**3) * (omega**3)
                    + 3 * I * zeta * omega * sympy.tan(theta / 2)
                    + 6
                ),
                (1 / 6)
                * (
                    (zeta**3) * (omega**3)
                    - 3 * zeta * omega * sympy.tan(theta / 2)
                    - 6
                ),
                0,
            ],
        ]
    )

    M3 = sympy.Matrix(
        [
            [I * sympy.sin(theta), -sympy.sin(theta), I * sympy.cos(theta) / omega],
            [
                -0.5
                * I
                * (
                    (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                    + 2 * sympy.cos(theta)
                ),
                -0.5 * (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                + sympy.cos(theta),
                -I * sympy.sin(theta) / omega,
            ],
            [
                -(1 / 6)
                * I
                * (
                    (zeta**3) * (omega**3)
                    + 3 * zeta * omega * sympy.tan(theta / 2)
                    + 6 * I
                ),
                (1 / 6)
                * (
                    -(zeta**3) * (omega**3)
                    + 3 * zeta * omega * sympy.tan(theta / 2)
                    - 6
                ),
                0,
            ],
        ]
    )
    M4 = sympy.Matrix(
        [
            [
                -sympy.sin(theta),
                -sympy.sin(theta),
                -sympy.cos(theta),
            ],
            [
                -sympy.cos(theta),
                -sympy.cos(theta),
                sympy.sin(theta),
            ],
            [
                I,
                1,
                0,
            ],
        ]
    )
    M5 = sympy.Matrix(
        [
            [
                -I * zeta * omega * (sympy.sin(theta / 2) ** 2),
                -zeta * omega * (sympy.sin(theta / 2) ** 2),
                -1,
            ],
            [
                -1 - 0.5 * I * zeta * theta * sympy.sin(theta),
                -1 - 0.5 * zeta * theta * sympy.sin(theta),
                0,
            ],
            [
                -I,
                -1,
                0,
            ],
        ]
    )
    M6 = sympy.Matrix(
        [
            [
                I * zeta * omega * (sympy.sin(theta / 2) ** 2),
                zeta * omega * (sympy.sin(theta / 2) ** 2),
                -1,
            ],
            [
                -1 + 0.5 * I * zeta * theta * sympy.sin(theta),
                -1 + 0.5 * zeta * theta * sympy.sin(theta),
                0,
            ],
            [
                I,
                1,
                0,
            ],
        ]
    )

    return [M1, M2, M3, M4, M5, M6]


if __name__ == "__main__":
    obj = get_equation_cross()
    filehandler = open("equations/cross_section", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    print("Cross Section equations saved")

    obj = get_equation_free_end()
    filehandler = open("equations/free_end", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    print("Free end equations saved")

    obj = get_equation_fixed_end()
    filehandler = open("equations/fixed_end", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    print("Fixed End equations saved")
