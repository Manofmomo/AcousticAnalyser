from sympy import symbols, Matrix, I, E
import logging
from typing import List
import numpy as np
import bisect 

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

        self.C = (self.youngs_modulus / self.density) ** 0.5
        self.K = (self.inertia / self.cross_section_area) ** 0.5

        self.constraint_count = 0
        self.constraint_ids = []

        self.force_counter=0
        self.force_pos = []

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

    def create_force_parameters(self) -> List[symbols]:

        g1_b_plus, g1_e_plus, g1_b_minus, g1_e_minus, g1_l_plus, g1_l_minus = symbols(
            "g_{f}b^+{i}, g_{f}e^+{i}, g_{f}b^-{i}, g_{f}e^-{i}, g_{f}l^+{i}, g_{f}l^-{i}".format(
                i=self.id, f = self.force_counter*2-1
            )
        )
        g1_plus = Matrix([g1_b_plus, g1_e_plus, g1_l_plus])
        g1_minus = Matrix([g1_b_minus, g1_e_minus, g1_l_minus])

        g2_b_plus, g2_e_plus, g2_b_minus, g2_e_minus, g2_l_plus, g2_l_minus = symbols(
            "g_{f}b^+{i}, g_{f}e^+{i}, g_{f}b^-{i}, g_{f}e^-{i}, g_{f}l^+{i}, g_{f}l^-{i}".format(
                i=self.id, f = self.force_counter*2
            )
        )
        g2_plus = Matrix([g2_b_plus, g2_e_plus, g2_l_plus])
        g2_minus = Matrix([g2_b_minus, g2_e_minus, g2_l_minus])
        
        return [g1_plus, g1_minus, g2_plus, g2_minus]


    def add_force(self, constraint_id: int, force: float, x: float):
        """Allows adding forced to members for forced vibration"""
        if x>=self.length or x<=0:
            raise ValueError("force lies outside the member")
        if x in self.force_pos:
            raise ValueError("point already has forces acting on it")
        
        self.force_counter = self.force_counter + 1
        # Correct for the sign convention
        if constraint_id == max(self.constraint_ids):
            x=self.length-x

        force_params = self.create_force_parameters()

        index = bisect.bisect_left(self.force_pos,x)
        self.force_pos.insert(index,x)
        self.params = self.params[:index*4+2] + force_params + self.params[index*4+2:]

    def get_propagation_matrix(self, w: float, lengths: List[float]):
        length = np.array(lengths)

        alpha = (w * self.K / self.C) ** 0.5
        beta = w * self.K / self.C
        L_bar = length / self.K

        propagation_matrix_finder = lambda x: np.array(
            [
                [np.e ** (-1j * alpha * x), 0, 0],
                [0, np.e ** (-alpha * x), 0],
                [0, 0, np.e ** (-1j * beta * x)],
            ]
        )
        propagation_matrix = np.stack([propagation_matrix_finder(x) for x in L_bar])
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

        self.a_plus = Matrix([a_b_plus, a_e_plus, a_l_plus])
        self.a_minus = Matrix([a_b_minus, a_e_minus, a_l_minus])

        self.b_plus: Matrix = Matrix([b_b_plus, b_e_plus, b_l_plus])
        self.b_minus: Matrix = Matrix([b_b_minus, b_e_minus, b_l_minus])
    
        self.params = [self.a_plus, self.a_minus, self.b_plus, self.b_minus]
        logger.debug(f"Parameters set for member {self.id}")

    def get_all_parameters(self) -> List[symbols]:
        # Flatten out parameters
        return [param for sublist in self.params for param in sublist]

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

        if len(self.force_pos)==0:
            propagation_matrix_subs = self.get_propagation_matrix(
                w=w, lengths=[self.length]
            )[0]
            matrix_forward = propagation_matrix_subs * self.a_plus - self.b_plus
            matrix_backward = propagation_matrix_subs * self.b_minus - self.a_minus
            logger.debug(f"Propagation for id:{self.id} calculated")
            eqns = matrix_forward.col_join(matrix_backward)
            return eqns
        
        else:
            for i in range(self.force_counter):
                x = self.force_pos[i]
                propagation_matrix_subs = self.get_propagation_matrix(
                w=w, lengths=[x]
                )[0]
                first_plus,first_minus,second_plus,second_minus, third_plus, third_minus = self.params[i*4:i*4+6]
                matrix_forward = propagation_matrix_subs * first_plus - second_plus
                matrix_backward = propagation_matrix_subs * second_minus - first_minus
                # matrix_force = 


    def get_non_dimensional_freq(self, w: float) -> float:
        omega = (w * self.K / self.C) ** 0.5
        return omega

    def get_deformation(
        self, w: float, lengths: List[float], subs_dict: dict, id: int
    ) -> List[np.ndarray]:
        lengths = np.array(lengths)

        a_minus, a_plus = self.get_parameters(w=w, id=id)
        # a_minus, a_plus = (self.a_minus, self.a_plus)

        propagation_matrix_subs = self.get_propagation_matrix(w=w, lengths=lengths)
        propagation_matrix_inv_subs = np.linalg.inv(propagation_matrix_subs)
        a_plus_subs = np.matrix(a_plus.subs(subs_dict), dtype=complex)
        a_minus_subs = np.matrix(a_minus.subs(subs_dict), dtype=complex)
        v = (
            np.matmul(np.array([1, 1, 0]), propagation_matrix_subs) * a_plus_subs
            + np.matmul(np.array([1, 1, 0]), propagation_matrix_inv_subs) * a_minus_subs
        )
        u = (
            np.matmul(np.array([0, 0, 1]), propagation_matrix_subs) * a_plus_subs
            + np.matmul(np.array([0, 0, 1]), propagation_matrix_inv_subs) * a_minus_subs
        )

        return [v, u]
