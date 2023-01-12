from sympy import symbols, Matrix, I, E
import logging
from typing import List

logger = logging.getLogger("acoustic_analyser")


class member:
    """This class emulates a member with all its physical properties"""

    def __init__(
        self,
        length: float,
        density: float,
        youngs_modulus: float,
        height: float,
        omega: symbols,
        id: int,
    ) -> None:
        if length <= 0:
            raise ValueError("length must be greater than 0")
        if density <= 0:
            raise ValueError("density must be greater than 0")
        if youngs_modulus <= 0:
            raise ValueError("youngs_modulus must be greater than 0")
        if height <= 0:
            raise ValueError("height must be greater than 0")

        self.length = length
        self.density = density
        self.youngs_modulus = youngs_modulus
        self.height = height
        self.id = id
        self.omega = omega

        self.cross_section_area = height * height
        self.inertia = (height**4) / 12

        self.constraint_count = 0
        self.constraint_ids = []

        self.set_parameters()

    def check_constraint_count(self) -> bool:
        """Each member can have only 2 constraints added to it at present"""
        if self.constraint_count < 2:
            return True
        else:
            return False

    def increment_constraint_count(self):
        self.constraint_count = self.constraint_count + 1

    def add_constraint(self, id):
        self.constraint_ids.append(id)

    def get_propagation_matrix(self, w: float, length: float):
        C = (self.youngs_modulus / self.density) ** 0.5
        K = (self.inertia / self.cross_section_area) ** 0.5

        alpha = (w * K / C) ** 0.5
        beta = w * K / C
        L_bar = length / K

        propagation_matrix = Matrix(
            [
                [E ** (-I * alpha * L_bar), 0, 0],
                [0, E ** (-alpha * L_bar), 0],
                [0, 0, E ** (-I * beta * L_bar)],
            ]
        )
        logger.debug(f"Propagation Matrix Calculated for member {self.id}")
        return propagation_matrix

    def set_parameters(self) -> None:
        a_b_plus, a_e_plus, a_b_minus, a_e_minus, a_l_plus, a_l_minus = symbols(
            "a_b^+{i}, a_e^+{i}, a_b^-{i}, a_e^-{i}, a_l^+{i}, a_l^-{i}".format(
                i=self.id
            )
        )
        b_b_plus, b_e_plus, b_b_minus, b_e_minus, b_l_plus, b_l_minus = symbols(
            "b_b^+{i}, b_e^+{i}, b_b^-{i}, b_e^-{i}, b_l^+{i}, b_l^-{i}".format(
                i=self.id
            )
        )
        self.params = [
            a_b_plus,
            a_e_plus,
            a_b_minus,
            a_e_minus,
            a_l_plus,
            a_l_minus,
            b_b_plus,
            b_e_plus,
            b_b_minus,
            b_e_minus,
            b_l_plus,
            b_l_minus,
        ]

        self.a_plus = Matrix([a_b_plus, a_e_plus, a_l_plus])
        self.a_minus = Matrix([a_b_minus, a_e_minus, a_l_minus])

        self.b_plus: Matrix = Matrix([b_b_plus, b_e_plus, b_l_plus])
        self.b_minus: Matrix = Matrix([b_b_minus, b_e_minus, b_l_minus])
        logger.debug(f"Parameters set for member {self.id}")

    def get_all_parameters(self) -> List[symbols]:
        return self.params

    def get_parameters(self, w: float, id: int) -> list:
        """This function gives back the set of parameters to be used.
        It corrects for the sign convention of the constraint when returning parameters
        Positive direction is from lower to higher constraint id
        """

        if id == max(self.constraint_ids):
            return [self.b_plus, self.b_minus]
        else:
            return [self.a_minus, self.a_plus]

    def get_equations(self, w: float):
        self.propagation_matrix_subs = self.get_propagation_matrix(w=w, length=self.length)
        matrix_forward = self.propagation_matrix_subs * self.a_plus - self.b_plus
        matrix_backward = self.propagation_matrix_subs * self.b_minus - self.a_minus
        logger.debug(f"Propagation for id:{self.id} calculated")
        eqns = matrix_forward.col_join(matrix_backward)
        return eqns

    def get_non_dimensional_freq(self, w: float) -> float:
        C = (self.youngs_modulus / self.density) ** 0.5
        K1 = (self.inertia / self.cross_section_area) ** 0.5
        omega = (w * K1 / C) ** 0.5
        return omega

    def get_deformation(self, w: float, length: float,subs_dict:dict, id: int) -> List[float]:
        if id == max(self.constraint_ids):
            length=self.length-length

        propagation_matrix_subs = self.get_propagation_matrix(w=w,length=length)
        propagation_matrix_inv_subs = propagation_matrix_subs.inv()
        a_plus_subs = self.a_plus.subs(subs_dict)
        a_minus_subs = self.a_minus.subs(subs_dict)

        v = (Matrix([1,1,0]).T*propagation_matrix_subs*a_plus_subs + Matrix([1,1,0]).T*propagation_matrix_inv_subs*a_minus_subs).evalf()
        u = (Matrix([0,0,1]).T*propagation_matrix_subs*a_plus_subs + Matrix([0,0,1]).T*propagation_matrix_inv_subs*a_minus_subs).evalf()

        return [v[0],u[0]]
