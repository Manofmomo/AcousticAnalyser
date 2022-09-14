from src.modules.member import member as member_type
from src.modules.rt_joint import get_rt_of_cross_section
from sympy import Matrix
from typing import List


class joint:
    force = 0.0
    k_spring = 0.0

    def __init__(self, members: list, id: int) -> None:
        for member in members:
            if type(member) != member_type:
                raise Exception("members must belong to class Member")
            if not member.check_constraint_count():
                raise Exception(
                    f"Member {member.member_number} already has 2 joints/BCs"
                )
        for member in members:
            member.increment_constraint_count()
            member.add_constraint(id=id)

        self.members: List[member_type] = members
        self.id = id

    def add_spring(self, k_spring: float) -> None:
        if k_spring <= 0:
            raise Exception("k_spring must be greater than 0")

        self.k_spring = k_spring

    def add_force(self, F: float) -> None:
        self.force = F


class two_member(joint):
    def __init__(
        self, theta: float, member_1: member_type, member_2: member_type, id: int
    ) -> None:
        self.theta = theta
        self.member_joint_type = []
        members = [member_1, member_2]
        super().__init__(members=members, id=id)
        self.reflection_matrix, self.transmission_matrix = get_rt_of_cross_section(
            m1=member_1, m2=member_2, theta=theta
        )

    def get_equations(self) -> list:
        # Assuming Wave equation to be A*e^(-I*k*x) + B*e^(-k*x) + C*e^(I*k*x) + D*e^(k*x)
        # Assuming
        c1, d1, a1, b1 = self.members[0].get_parameters(id=self.id)
        a2, b2, _, _ = self.members[1].get_parameters(id=self.id)

        matrix_reflect = self.reflection_matrix * Matrix([a1, b1]) - Matrix([c1, d1])
        matrix_transmit = self.transmission_matrix * Matrix([a1, b1]) - Matrix([a2, b2])
        eqns = matrix_reflect.col_join(matrix_transmit)

        return eqns
