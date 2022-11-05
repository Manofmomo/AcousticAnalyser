import sympy
from sympy import I
import pickle
from typing import Tuple
import logging

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
    """Gets the M0-M6 of a joint assuming Wave equation to be A*e^(-I*k*x) + B*e^(-k*x) + C*e^(I*k*x) + D*e^(k*x)"""

    w = sympy.symbols("w")
    density1, area1, E1, I1, H1 = sympy.symbols("density1, area1, E1, I1, H1")

    C = (E1 / density1) ** 0.5
    K1 = (I1 / area1) ** 0.5

    zeta = H1 / K1
    omega = (w * K1 / C) ** 0.5
    theta = sympy.symbols("theta")

    M1 = sympy.Matrix(
        [
            [0, 0, zeta * omega - I / omega],
            [zeta * omega - I, zeta * omega + 1, 0],
            [
                1 + 0.5 * I * zeta * omega * sympy.tan(theta / 2),
                -1 - 0.5 * I * zeta * omega * sympy.tan(theta / 2),
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
            [I * sympy.sin(theta), -sympy.sin(theta), -I * sympy.cos(theta) / omega],
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
                -1 - 0.5 * I * zeta * omega * sympy.sin(theta),
                -1 - 0.5 * zeta * omega * sympy.sin(theta),
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
                -1 + 0.5 * I * zeta * omega * sympy.sin(theta),
                -1 + 0.5 * zeta * omega * sympy.sin(theta),
                0,
            ],
            [
                I,
                1,
                0,
            ],
        ]
    )

    N1 = sympy.Matrix(
        [
            [I * sympy.sin(theta), -sympy.sin(theta), -I * sympy.cos(theta) / omega],
            [-I * sympy.cos(theta), sympy.cos(theta), -I * sympy.sin(theta) / omega],
            [
                1 - 0.5 * I * zeta * omega * sympy.tan(theta / 2),
                1 + 0.5 * 1 * zeta * omega * sympy.tan(theta / 2),
                0,
            ],
        ]
    )

    N2 = sympy.Matrix(
        [
            [0, 0, -zeta * omega + I / omega],
            [
                -0.5 * I * (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                - zeta * omega
                + I,
                -0.5 * 1 * (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                - zeta * omega
                - 1,
                0,
            ],
            [
                -(1 / 6)
                * I
                * (
                    (zeta * omega) ** 3
                    + 3 * zeta * omega * sympy.tan(theta / 2)
                    - 6 * I
                ),
                (1 / 6)
                * 1
                * (
                    -((zeta * omega) ** 3) + 3 * zeta * omega * sympy.tan(theta / 2) + 6
                ),
                0,
            ],
        ]
    )

    N3 = sympy.Matrix(
        [
            [0, 0, -zeta * omega - I / omega],
            [
                0.5
                * I
                * (
                    (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                    + 2 * I * zeta * omega
                    - 2
                ),
                0.5 * 1 * (zeta**2) * (omega**2) * sympy.tan(theta / 2)
                - zeta * omega
                + 1,
                0,
            ],
            [
                (1 / 6)
                * I
                * (
                    (zeta * omega) ** 3
                    + 3 * zeta * omega * sympy.tan(theta / 2)
                    + 6 * I
                ),
                (1 / 6)
                * 1
                * ((zeta * omega) ** 3 - 3 * zeta * omega * sympy.tan(theta / 2) + 6),
                0,
            ],
        ]
    )

    N4 = sympy.Matrix(
        [
            [0, 0, 1],
            [1, 1, 0],
            [-I, -1, 0],
        ]
    )

    N5 = sympy.Matrix(
        [
            [
                sympy.sin(theta) + I * zeta * omega * (sympy.sin(theta / 2)) ** 2,
                0.5 * sympy.sin(theta) * (zeta * omega * sympy.tan(theta / 2) + 2),
                sympy.cos(theta),
            ],
            [
                sympy.cos(theta) + 0.5 * I * zeta * omega * sympy.sin(theta),
                0.5 * zeta * omega * sympy.sin(theta) + sympy.cos(theta),
                -sympy.sin(theta),
            ],
            [I, 1, 0],
        ]
    )

    N6 = sympy.Matrix(
        [
            [
                sympy.sin(theta) - I * zeta * omega * (sympy.sin(theta / 2)) ** 2,
                sympy.sin(theta) - zeta * omega * (sympy.sin(theta / 2)) ** 2,
                sympy.cos(theta),
            ],
            [
                sympy.cos(theta) - 0.5 * I * zeta * omega * sympy.sin(theta),
                sympy.cos(theta) - 0.5 * 1 * zeta * omega * sympy.sin(theta),
                -sympy.sin(theta),
            ],
            [-I, -1, 0],
        ]
    )

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
