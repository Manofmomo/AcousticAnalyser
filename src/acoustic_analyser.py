from itertools import count
from src.modules import joint, bc
from sympy.core.add import Add as equation_type
from sympy import symbols, Matrix
from src.modules.member import member as member_type
from json import load as json_load
from csv import reader as csv_load


def get_coefficients(eqn: equation_type, params: list) -> list:
    # This is intended to extract the coefficients given an equation and the variables
    eqn = eqn.expand()
    coeffs = []
    for param in params:
        coeffs.append(eqn.coeff(param))
    return coeffs


class frame:
    # This class emulates a mechanical frame consisting of joints and members
    # It is the class the user directly interacts with
    def __init__(self) -> None:
        self.members = {}
        self.constraints = []
        self.constraints_count = -1
        self.omega = symbols("w")

    @classmethod
    def from_file(cls, member_file: str, constraint_file: str):
        obj = cls()

        with open(member_file, "r") as jsonfile:
            member_dict = json_load(jsonfile)
            counter = 0
            for member_id, member_deets in member_dict.items():
                assert int(member_id) == counter, "ID format is not followed"
                obj.add_member(id=int(member_id), **member_deets)
                counter = counter + 1

        with open(constraint_file, "r") as csvfile:
            constraints = csv_load(csvfile)
            for m1_id, row in enumerate(constraints):
                for m2_id, val in enumerate(row):
                    if int(val) == -1:
                        continue

                    print(m1_id, m2_id)
                    # For >2 member joints, we can add a condition to check how many joints are there and appropriately select function
                    obj.two_member_joint(
                        theta=float(val), member_1_id=int(m1_id), member_2_id=int(m2_id)
                    )

        return obj

    def add_member(
        self,
        length: float,
        density: float,
        cross_section_area: float,
        youngs_modulus: float,
        inertia: float,
        id: int,
    ) -> member_type:
        # This function allows the user to add a member to the structure
        member_obj = member_type(
            length=length,
            density=density,
            cross_section_area=cross_section_area,
            youngs_modulus=youngs_modulus,
            inertia=inertia,
            omega=self.omega,
            id=id,
        )
        self.members[id] = member_obj
        return member_obj

    def _add_constraint(func):
        # This is a decorator function defined to ease the process of development
        # It handles the adding of a constraint object to the structure
        def inner1(self, *args, **kwargs):
            constraint_obj = func(self, *args, **kwargs)
            self.constraints.append(constraint_obj)
            self.constraints_count = self.constraints_count + 1
            return constraint_obj

        return inner1

    def _get_constraint_id(self):
        # Gets the constraint ID that would be used if the checks succeed
        return self.constraints_count + 1

    # Here the joints/BCs are defined
    @_add_constraint
    def free_end(self, member_id: int) -> bc:
        free_end_obj = bc.free_end(
            member=self.members[member_id], id=self._get_constraint_id()
        )
        return free_end_obj

    @_add_constraint
    def fixed_end(self, member_id: int) -> bc:
        fixed_end_obj = bc.fixed_end(
            member=self.members[member_id], id=self._get_constraint_id()
        )
        return fixed_end_obj

    @_add_constraint
    def two_member_joint(
        self, theta: float, member_1_id: int, member_2_id: int
    ) -> joint:
        rigid_joint_obj = joint.two_member(
            theta=theta,
            member_1=self.members[member_1_id],
            member_2=self.members[member_2_id],
            id=self._get_constraint_id(),
        )
        return rigid_joint_obj

    def get_equation_matrix(self):
        # This function is responsible for collecting the equations from the constraints and constructing the desired matrix from them
        eqns = []
        for constraint in self.constraints:
            eqns.extend(constraint.get_equations())
        return eqns
