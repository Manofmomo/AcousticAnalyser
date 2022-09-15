from src.modules import joint, bc
from sympy.core.add import Add as equation_type
from sympy import symbols, Matrix
from src.modules.member import member as member_type


def get_coefficients(eqn: equation_type, params: list) -> list:
    # This is intended to extract the coefficients given an equation and the variables
    eqn = eqn.expand()
    coeffs = []
    for param in params:
        coeffs.append(eqn.coeff(param))
    return coeffs


class structure:
    # This class emulates a mechanical structure consisting of joints and members
    # It is the class the user directly interacts with
    def __init__(self) -> None:
        self.members = []
        self.member_count = -1
        self.constraints = []
        self.constraints_count = -1
        self.omega = symbols("w")

    def add_member(
        self,
        length: float,
        rho: float,
        cross_section_area: float,
        youngs_modulus: float,
        inertia: float,
    ) -> member_type:
        # This function allows the user to add a member to the structure
        member_obj = member_type(
            length=length,
            rho=rho,
            cross_section_area=cross_section_area,
            youngs_modulus=youngs_modulus,
            inertia=inertia,
            omega=self.omega,
            id=self.member_count,
        )
        self.members.append(member_obj)
        self.member_count = self.member_count + 1
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
    def free_end(self, member: member_type) -> bc:
        free_end_obj = bc.free_end(member=member, id=self._get_constraint_id())
        return free_end_obj

    @_add_constraint
    def fixed_end(self, member: member_type) -> bc:
        fixed_end_obj = bc.fixed_end(member=member, id=self._get_constraint_id())
        return fixed_end_obj

    @_add_constraint
    def two_member_joint(
        self, theta: float, member_1: member_type, member_2: member_type
    ) -> joint:
        rigid_joint_obj = joint.two_member(
            theta=theta,
            member_1=member_1,
            member_2=member_2,
            id=self._get_constraint_id(),
        )
        return rigid_joint_obj

    def get_equation_matrix(self):
        # This function is responsible for collecting the equations from the constraints and constructing the desired matrix from them
        eqns = []
        params = []
        lhs_list = []
        rhs_list = []
        for constraint in self.constraints:
            eqns.extend(constraint.get_equations())

        # for member in self.members:
        #     params.extend(member.get_parameters())

        # for lhs in eqns:
        #     lhs_list.append(get_coefficients(lhs, params))
        #     # rhs_list.append(rhs)

        return eqns
