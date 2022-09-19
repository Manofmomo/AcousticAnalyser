import sympy
from sympy import I
import pickle

# This initial code is developed using Euler-Bernoulli beam mode


def get_equation_free_end():

    x = sympy.symbols("x")
    w = sympy.symbols("w")
    rho1, Area1, E1, I1 = sympy.symbols("rho1, Area1, E1, I1")

    k1 = (rho1 * Area1 / (E1 * I1)) ** (0.25) * (w**0.5)

    c1, d1 = sympy.symbols("c_traveling, d_traveling")
    c2, d2 = sympy.symbols("c_evanescent, d_evanescent")
    params = [(c1, d1), (c2, d2)]

    travelling_wave = sympy.exp(-I * k1 * x)
    evanescent_wave = sympy.exp(-1 * k1 * x)
    waves = [travelling_wave, evanescent_wave]

    final_equations = []

    for i in range(2):
        incoming_wave = waves[i]
        c, d = params[i]

        v = incoming_wave + c * sympy.exp(I * k1 * x) + d * sympy.exp(k1 * x)
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
    rho1, Area1, E1, I1 = sympy.symbols("rho1, Area1, E1, I1")

    k1 = (rho1 * Area1 / (E1 * I1)) ** (0.25) * (w**0.5)

    c1, d1 = sympy.symbols("c_traveling, d_traveling")
    c2, d2 = sympy.symbols("c_evanescent, d_evanescent")
    params = [(c1, d1), (c2, d2)]

    travelling_wave = sympy.exp(-I * k1 * x)
    evanescent_wave = sympy.exp(-1 * k1 * x)
    waves = [travelling_wave, evanescent_wave]

    final_equations = []

    for i in range(2):
        incoming_wave = waves[i]
        c, d = params[i]

        v = incoming_wave + c * sympy.exp(I * k1 * x) + d * sympy.exp(k1 * x)
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

    x = sympy.symbols("x")
    w = sympy.symbols("w")
    rho1, Area1, E1, I1 = sympy.symbols("rho1, Area1, E1, I1")
    rho2, Area2, E2, I2 = sympy.symbols("rho2, Area2, E2, I2")

    k1 = (rho1 * Area1 / (E1 * I1)) ** (0.25) * (w**0.5)
    k2 = (rho2 * Area2 / (E2 * I2)) ** (0.25) * (w**0.5)

    a1, b1, c1, d1 = sympy.symbols(
        "a_traveling,  b_traveling, c_traveling, d_traveling"
    )
    a2, b2, c2, d2 = sympy.symbols(
        "a_evanescent, b_evanescent, c_evanescent, d_evanescent"
    )
    params = [(a1, b1, c1, d1), (a2, b2, c2, d2)]

    travelling_wave = sympy.exp(-I * k1 * x)
    evanescent_wave = sympy.exp(-1 * k1 * x)
    waves = [travelling_wave, evanescent_wave]

    final_equations = []

    for i in range(2):
        incoming_wave = waves[i]
        a, b, c, d = params[i]

        v1 = incoming_wave + c * sympy.exp(I * k1 * x) + d * sympy.exp(k1 * x)
        v2 = a * sympy.exp(-I * k2 * x) + b * sympy.exp(-k2 * x)
        v = sympy.Piecewise((v1, x < 0), (v2, x > 0))
        vdiffs = [None] * 4
        vdiffs[0] = v
        for i in range(3):
            vdiffs[i + 1] = sympy.diff(vdiffs[i], x)

        eq1 = sympy.limit(vdiffs[0], x, 0, dir="-") - sympy.limit(
            vdiffs[0], x, 0, dir="+"
        )
        eq2 = sympy.limit(vdiffs[1], x, 0, dir="-") - sympy.limit(
            vdiffs[1], x, 0, dir="+"
        )
        eq3 = E1 * I1 * sympy.limit(vdiffs[2], x, 0, dir="-") - E2 * I2 * sympy.limit(
            vdiffs[2], x, 0, dir="+"
        )
        eq4 = E1 * I1 * sympy.limit(vdiffs[3], x, 0, dir="-") - E2 * I2 * sympy.limit(
            vdiffs[3], x, 0, dir="+"
        )

        final_equations.append([eq1, eq2, eq3, eq4])

    return final_equations


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
