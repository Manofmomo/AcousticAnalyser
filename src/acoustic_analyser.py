from src.modules import joint, bc
from src.modules.bc import bc as bc_type
from src.modules.joint import joint as joint_type
from sympy.core.add import Add as equation_type
from sympy import symbols, Matrix
from src.modules.member import member as member_type
from json import load as json_load
from csv import reader as csv_load
import logging
from typing import Dict, List, Union
import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse.linalg import eigs
import itertools
from scipy.optimize import newton

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s:%(message)s"
)
logger = logging.getLogger("acoustic_analyser")


def get_coefficient_matrix(eqns: equation_type, params: list) -> Matrix:
    """This is intended to extract the coefficients from eqns and build a matrix"""
    coeff_list = []
    for eqn in eqns:
        eqn = eqn.expand()
        coeffs = []
        for param in params:
            coeffs.append(complex(eqn.coeff(param)))
        coeff_list.append(coeffs)
    return np.array(coeff_list)


class frame:
    """This class emulates a mechanical frame consisting of joints and members
    It is the class the user directly interacts with
    """

    def __init__(self, debug: bool = False) -> None:
        if debug:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        self.debug = debug
        self.members: Dict[int, member_type] = {}
        self.constraints: Dict[int, Union[bc_type, joint_type]] = {}
        self.constraints_count = -1
        self.omega = symbols("w")

    @classmethod
    def from_file(cls, member_file: str, constraint_file: str, debug: bool = False):
        obj = cls(debug)

        with open(member_file, "r") as jsonfile:
            member_dict = json_load(jsonfile)
            counter = 0
            for member_id, member_deets in member_dict.items():
                assert int(member_id) == counter, "ID format is not followed"
                obj.add_member(id=int(member_id), **member_deets)
                counter = counter + 1
        logger.debug("Created all Members")
        with open(constraint_file, "r") as csvfile:
            constraints = csv_load(csvfile)
            for m1_id, row in enumerate(constraints):
                for m2_id, val in enumerate(row):
                    if int(val) == -1:
                        continue
                    # For >2 member joints, we can add a condition to check how many joints are there and appropriately select function
                    obj.two_member_joint(
                        theta=float(val), member_1_id=int(m1_id), member_2_id=int(m2_id)
                    )
        logger.debug("Added all Constraints")

        return obj

    def add_member(
        self,
        length: float,
        density: float,
        youngs_modulus: float,
        height: float,
        id: int,
    ) -> member_type:
        """This function allows the user to add a member to the structure"""
        member_obj = member_type(
            length=length,
            density=density,
            youngs_modulus=youngs_modulus,
            height=height,
            omega=self.omega,
            id=id,
        )
        self.members[id] = member_obj
        logger.debug(f"Member added with id {id}")
        return member_obj

    def _add_constraint(func):
        """This is a decorator function defined to ease the process of development
        It handles the adding of a constraint object to the structure
        """

        def inner1(self, *args, **kwargs):
            constraint_obj = func(self, *args, **kwargs)
            self.constraints[self._get_constraint_id()] = constraint_obj
            self.constraints_count = self.constraints_count + 1
            return constraint_obj

        return inner1

    def _get_constraints(self):
        return self.constraints

    def _get_constraint_id(self):
        """Gets the constraint ID that would be used if the checks succeed"""
        return self.constraints_count + 1

    # Here the joints/BCs are defined
    @_add_constraint
    def free_end(self, member_id: int) -> bc:
        id = self._get_constraint_id()
        free_end_obj = bc.free_end(member=self.members[member_id], id=id)
        logger.debug(f"Free End added to member {member_id} with id {id}")
        return free_end_obj

    @_add_constraint
    def fixed_end(self, member_id: int) -> bc:
        id = self._get_constraint_id()
        fixed_end_obj = bc.fixed_end(member=self.members[member_id], id=id)
        logger.debug(f"Fixed End added to member {member_id} with id {id}")
        return fixed_end_obj

    @_add_constraint
    def two_member_joint(
        self, theta: float, member_1_id: int, member_2_id: int
    ) -> joint:
        id = self._get_constraint_id()
        rigid_joint_obj = joint.two_member(
            theta=theta,
            member_1=self.members[member_1_id],
            member_2=self.members[member_2_id],
            id=id,
        )
        logger.debug(
            f"Two member joint added b/w {member_1_id} and {member_2_id} with id {id}"
        )
        return rigid_joint_obj

    def _set_params(self) -> None:
        self.params = []
        for member in self.members.values():
            self.params.extend(member.get_all_parameters())

    def get_equation_matrix(self, w: float):
        """This function is responsible for collecting the equations from the constraints and constructing the desired matrix from them"""
        self._set_params()

        eqns = Matrix([])
        for constraint in self.constraints.values():
            eqns = eqns.col_join(constraint.get_equations(w=w))
        for member in self.members.values():
            eqns = eqns.col_join(member.get_equations(w=w))
        logger.debug("All Equations Fetched")
        self.eqn_matrix = eqns
        for eqn in eqns:
            logger.debug(eqn.expand())
        coeff_matrix = get_coefficient_matrix(eqns=eqns, params=self.params)
        logger.debug("Equation Coefficent Matrix generated")
        return coeff_matrix

    def get_determinant(self, w: float) -> complex:
        """Returns the determinant of the A matrix given omega"""
        det = np.linalg.det(self.get_equation_matrix(w=w))
        print(f"Determinant: {det}")
        return det

    def get_frequency_graph(
        self, lower_limit: float, upper_limit: float, step_size: float
    ) -> np.array:
        "Creates a graph of the real and imaginary components of the determinant"
        freq = np.arange(lower_limit, upper_limit, step_size)
        output = np.zeros(len(freq), dtype=np.complex_)
        omega = freq * np.pi * 2

        for i in range(len(freq)):
            output[i] = self.get_determinant(w=omega[i])
            print(f"{i}/{len(freq)}\t", end="\r")

        plt.plot(freq, np.real(output))
        plt.plot(freq, np.imag(output))
        plt.plot([lower_limit, upper_limit], [0, 0], "--")
        plt.legend(["Real", "Imaginary"])
        plt.show()
        plt.figure()
        plt.plot(freq, np.abs(output))
        plt.legend(["Abs"])
        plt.show()
        return output

    def get_natural_frequency_newton(
        self,
        initial_guess: float,
        tol: float = 1e-09,
        max_iter: int = 100,
    ) -> float:
        """Takes in the a initial guess and uses newton raphson method to solve for natural frequency"""
        get_determinant = lambda x: self.get_determinant(w=np.abs(x) * 2 * np.pi)
        return np.abs(newton(get_determinant, initial_guess, tol=tol, maxiter=max_iter))

    def get_natural_frequency_bisect(
        self,
        lower_limit: float,
        upper_limit: float,
        precision: float = 0.01,
        max_iter: int = 100,
    ) -> float:
        """Takes in the upper and lower limit of frequency between which it searches for the natural frequency of the structure.
        The natural frequency must lie between the two limits.
        """
        freq_1 = lower_limit
        freq_2 = upper_limit

        iter = 0
        Handle = True
        omega = lambda freq: freq * np.pi * 2

        def sign_check(output_1: complex, output_2: complex) -> bool:
            real_sign = np.sign(np.real(output_1) / np.real(output_2))
            imag_sign = np.sign(np.imag(output_1) / np.imag(output_2))
            if real_sign <= 0 and imag_sign < 0:
                return True
            else:
                return False

        output_1 = self.get_determinant(w=omega(freq=freq_1))
        output_2 = self.get_determinant(w=omega(freq=freq_2))

        while Handle:
            iter = iter + 1
            if iter > max_iter:
                Handle = False

            freq_12 = (freq_1 + freq_2) / 2
            output_12 = self.get_determinant(w=omega(freq=freq_12))

            if np.isclose(abs(output_12), 0) or (freq_2 - freq_1) < precision:
                Handle = False

            if sign_check(output_1, output_12):
                freq_2 = freq_12
            elif sign_check(output_2, output_12):
                freq_1 = freq_12
            else:
                logger.error("There is no solution between these points")
                Handle = False
                return None

            output_1 = self.get_determinant(w=omega(freq=freq_1))
            output_2 = self.get_determinant(w=omega(freq=freq_2))

            print(
                f"Iteration: {iter}, Det: {abs(output_12)}, freq: {freq_12}\t", end="\r"
            )

        return (freq_1 + freq_2) / 2

    def get_params_solution(self, natural_freq, atol=1e-05, rtol=1e-05) -> Dict:
        """Using the natural frequncy and computes the parameters by finding eigenvector for the eigenvalue 0"""
        # matrix = self.get_equation_matrix(w=2 * np.pi * natural_freq)
        # val, vec = eigs(matrix, k=1, sigma=0)
        # val = val[0]

        val, vec = np.linalg.eig(self.get_equation_matrix(w=2 * np.pi * 13.742))
        solns = vec[:, np.abs(val).argmin()]

        assert (
            np.isclose(val[np.abs(val).argmin()], 0, atol=atol, rtol=rtol) == True
        ), "The value provided is not a natural frequency"

        self.params_subs = dict(zip(self.params, solns / 2))
        return self.params_subs
