from acoustic_analyser.modules import joint, bc
from acoustic_analyser.modules.bc import bc as bc_type
from acoustic_analyser.modules.bc import fixed_end as fixed_end_type
from acoustic_analyser.modules.joint import joint as joint_type
from sympy.core.add import Add as equation_type
from sympy import symbols, Matrix
from acoustic_analyser.modules.member import member as member_type
from json import load as json_load
from csv import reader as csv_load
import logging
from typing import Dict, List, Union
import numpy as np
from matplotlib import pyplot as plt
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
        """Defines a frame from a constraint file and a member file."""
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
        cross_section_area: float,
        inertia: float,
        id: int,
    ) -> member_type:
        """This function allows the user to add a member to the structure"""
        member_obj = member_type(
            length=length,
            density=density,
            youngs_modulus=youngs_modulus,
            cross_section_area=cross_section_area,
            inertia=inertia,
            omega=self.omega,
            height=height,
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

    def get_determinant(self, w: float, print_det: bool = False) -> complex:
        """Returns the determinant of the A matrix given omega"""
        det = np.linalg.det(self.get_equation_matrix(w=w))
        if print_det:
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
            print(f"{i}/{len(freq)}", " " * 10, end="\r")

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
        print_det: bool = True,
    ) -> float:
        """Takes in the a initial guess and uses newton raphson method to solve for natural frequency"""
        get_determinant = lambda x: self.get_determinant(
            w=np.abs(x) * 2 * np.pi, print_det=print_det
        )
        return np.abs(newton(get_determinant, initial_guess, tol=tol, maxiter=max_iter))

    def get_natural_frequency(
        self,
        n: int = 1,
        lower_limit: float = 0,
        upper_limit: float = 1000,
        step_size: float = 1,
        atol=1e-08,
        rtol=1e-05,
    ) -> List[float]:
        """Returns the first n natural frequencies by using the bisect method at each step between the lower and upper limit"""
        natural_frequencies = []

        freq_1 = lower_limit
        freq_2 = lower_limit + step_size
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

        # Implemented the bisect method
        while n > 0 and freq_1 < upper_limit:
            print(f"Checking region [{freq_1}, {freq_2}]", " " * 10, end="\r")
            if sign_check(output_1, output_2):
                print(f"\nPossible root in region [{freq_1}, {freq_2}]")
                max_iter = 100
                root_found = False
                counter = 0
                x1 = freq_1
                x2 = freq_2
                while counter < max_iter:
                    counter = counter + 1
                    x = (x1 + x2) / 2
                    output = self.get_determinant(w=omega(freq=x))
                    if np.isclose(output, 0, atol=atol, rtol=rtol):
                        natural_frequencies.append(x)
                        n = n - 1
                        root_found = True
                        break

                    if sign_check(output_1, output):
                        x2 = x
                    elif sign_check(output_2, output):
                        x1 = x
                    else:
                        break
                    print(
                        f"Iteration: {counter}, Det: {abs(output)}, freq: {x}", end="\r"
                    )
                print(" " * 100, end="\r")
                if not root_found:
                    print(f"Natural frequency not found in range, continuing")
                else:
                    print(f"Natural frequency found at {x}")
            freq_1 = freq_2
            freq_2 = freq_2 + step_size
            output_1 = output_2
            output_2 = self.get_determinant(w=omega(freq=freq_2))

        return natural_frequencies

    def get_params_solution(self, natural_freq, atol=1e-05, rtol=1e-05) -> Dict:
        """Using the natural frequncy and computes the parameters by finding eigenvector for the eigenvalue 0"""

        val, vec = np.linalg.eig(self.get_equation_matrix(w=2 * np.pi * natural_freq))
        solns = vec[:, np.abs(val).argmin()]

        assert (
            np.isclose(val[np.abs(val).argmin()], 0, atol=atol, rtol=rtol) == True
        ), "The value provided is not a natural frequency"

        self.params_subs = dict(zip(self.params, solns))
        return self.params_subs

    def get_mode_shape(
        self,
        n: int = None,
        natural_freq: float = None,
        origin_constraint_id: int = None,
        step_size: float = 0.01,
        scaling_factor: float = 0.5,
    ) -> np.array:
        """
        Finds the mode shape of the frame and plots it.
        Origin constraint id is the id of the constraint (joint/bc) which the user would like at origin. In general, it can be left blank.
        If origin constraint id is not provided, it fixed a fixed_end as the origin.
        """
        if (n is None) and (natural_freq is None):
            logger.error(
                "Either n or natural frequency must be provided in order to find the mode shape"
            )
            return None
        if natural_freq is None:
            print("Beginning search for Natural Frequency")
            natural_freqs = self.get_natural_frequency(n=n)
            natural_freq = natural_freqs[-1]
            print(f"Natural Frequency for {n} mode shape is {natural_freq}")
        omega = natural_freq * np.pi * 2

        for constraint in self.constraints.values():
            if type(constraint) == fixed_end_type:
                origin_constraint_id = constraint.id
                break
        if origin_constraint_id == None:
            logger.error(
                "No Fixed End found in figure, please specify origin_constraint_id"
            )
            return None

        _ = self.get_params_solution(natural_freq=natural_freq)
        constraint_curr = self.constraints[origin_constraint_id]
        member_curr = constraint_curr.members[0]

        # Initialising values
        offset = np.array([0, 0])
        angle = 0
        Handle = True
        mode_shape_positive = [np.matrix([[0], [0]]).reshape(1, 2)]
        mode_shape_negative = [np.matrix([[0], [0]]).reshape(1, 2)]
        original_shape = [np.matrix([[0], [0]]).reshape(1, 2)]
        members_completed = set()
        constraints_completed = set()

        while Handle:
            # Set of completed members and constraints to prevent loop
            members_completed.add(member_curr.id)
            constraints_completed.add(constraint_curr.id)

            # Points along memmber in its local coordinate system
            x = np.append(
                np.arange(step_size, member_curr.length, step_size), member_curr.length
            )
            v, u = member_curr.get_deformation(
                w=omega, lengths=x, id=constraint_curr.id, subs_dict=self.params_subs
            )
            v = v * scaling_factor
            u = u * scaling_factor

            x_deformed_positive = np.real(u).reshape(len(x)) + x
            y_deformed_positive = np.real(v).T

            x_deformed_negative = -np.real(u).reshape(len(x)) + x
            y_deformed_negative = -np.real(v).T

            points_deformed_positive = np.stack(
                [x_deformed_positive, y_deformed_positive]
            ).T.reshape(len(x), 2)
            points_deformed_negative = np.stack(
                [x_deformed_negative, y_deformed_negative]
            ).T.reshape(len(x), 2)
            points_original = np.stack([x, np.zeros(len(x))]).T.reshape(len(x), 2)

            # Converting local coordinates to global coordinates
            rotation_matrix = np.array(
                [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
            )
            deformation_rotated_positive = np.matmul(
                rotation_matrix, points_deformed_positive.T
            ).T
            deformation_rotated_translated_positive = (
                deformation_rotated_positive + offset
            )

            deformation_rotated_negative = np.matmul(
                rotation_matrix, points_deformed_negative.T
            ).T
            deformation_rotated_translated_negative = (
                deformation_rotated_negative + offset
            )

            original_rotated = np.matmul(rotation_matrix, points_original.T).T
            original_rotated_translated = original_rotated + offset

            mode_shape_positive.append(deformation_rotated_translated_positive)
            mode_shape_negative.append(deformation_rotated_translated_negative)
            original_shape.append(original_rotated_translated)

            offset = offset + np.array(
                np.matmul(rotation_matrix, np.array([member_curr.length, 0])).T
            ).reshape(
                2,
            )

            constraint_ids = member_curr.constraint_ids
            if constraint_ids[0] not in constraints_completed:
                constraint_curr = self.constraints[constraint_ids[0]]
            else:
                constraint_curr = self.constraints[constraint_ids[1]]

            angle = angle + constraint_curr.theta
            angle = np.mod(angle, 2 * np.pi)

            members = (
                constraint_curr.members
            )  # Will have to change for a 3 member joint

            if members[0].id not in members_completed:
                member_curr = members[0]
            elif len(members) > 1 and members[1] not in members_completed:
                member_curr = members[1]
            else:
                Handle = False

        mode_shape_positive = np.concatenate(mode_shape_positive)
        mode_shape_negative = np.concatenate(mode_shape_negative)
        original_shape = np.concatenate(original_shape)

        plt.plot(mode_shape_positive[:, 0], mode_shape_positive[:, 1], "b")
        plt.plot(original_shape[:, 0], original_shape[:, 1], "r--")
        plt.plot(mode_shape_negative[:, 0], mode_shape_negative[:, 1], "b")
        plt.legend(["Mode Shape", "Original Shape"])
        plt.xticks([])
        plt.yticks([])

        return (
            np.array(mode_shape_positive),
            np.array(mode_shape_negative),
            np.array(original_shape),
        )
