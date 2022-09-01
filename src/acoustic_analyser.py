from src.modules import joint, bc
from sympy.core.add import Add as equation_type
from sympy import symbols, Matrix
from src.modules.member import member as member_type


def get_coefficients(eqn: equation_type, params: list) -> list:
    eqn = eqn.expand()
    coeffs = []
    for param in params:
        coeffs.append(eqn.coeff(param))
    return coeffs


class structure:
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
    ) -> tuple:
        member_obj = member_type(
            length=length,
            rho=rho,
            cross_section_area=cross_section_area,
            youngs_modulus=youngs_modulus,
            inertia=inertia,
            omega=self.omega,
            member_number=self.member_count,
        )
        self.members.append(member_obj)
        self.member_joint_count[member_obj] = 0
        self.member_count = self.member_count + 1
        return (member_obj, self.member_count)

    def _add_constraint(func):
        def inner1(self, *args, **kwargs):
            constraint_obj = func(*args, **kwargs)
            self.constraints.append(constraint_obj)
            self.constraints_count = self.constraints_count + 1
            return (constraint_obj, self.constraints_count)

        return inner1

    @_add_constraint
    def free_end(self, member: member_type) -> tuple:
        free_end_obj = bc.free_end(member=member)
        return free_end_obj

    @_add_constraint
    def fixed_end(self, member: member_type) -> tuple:
        fixed_end_obj = bc.fixed_end(member=member)
        return fixed_end_obj

    @_add_constraint
    def two_member_joint(
        self, theta: float, member_1: member_type, member_2: member_type
    ) -> tuple:
        rigid_joint_obj = joint.two_member(
            theta=theta, member_1=member_1, member_2=member_2
        )
        return rigid_joint_obj

    def get_equation_matrix(self):
        eqns = []
        params = []
        lhs_list = []
        rhs_list = []
        for joint in self.joints:
            eqns.extend(joint.get_equations())

        for member in self.members:
            params.extend(member.get_parameters())

        for lhs, rhs in eqns:
            lhs_list.append(self.get_coefficients(lhs, params))
            rhs_list.append(rhs)

        return (Matrix(lhs_list), Matrix(rhs_list))
