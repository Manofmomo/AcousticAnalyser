from src.modules.member import member as member_type
from src.modules.rt_bc import get_r_of_fixed_end, get_r_of_free_end
from sympy import Matrix


class bc:
    def __init__(self, member: member_type, id: int) -> None:
        if type(member) != member_type:
            raise Exception("members must belong to class Member")
        if not member.check_constraint_count():
            raise Exception(f"Member {member.member_number} already has 2 joints/BCs")
        member.increment_constraint_count()
        member.add_constraint(id=id)
        self.member = member
        self.id = id


class free_end(bc):
    def __init__(self, member: member_type, id: int) -> None:
        super().__init__(member=member, id=id)
        self.reflection_matrix = get_r_of_free_end(m1=member)

    def get_equations(self) -> list:  # TODO add type
        # Assuming Wave equation to be A*e^(-I*k*x) + B*e^(-k*x) + C*e^(I*k*x) + D*e^(k*x)
        a1, b1, c1, d1 = self.member.get_parameters(id=self.id)

        eqns = self.reflection_matrix * Matrix([c1, d1]) - Matrix([a1, b1])

        return eqns


class fixed_end(bc):
    def __init__(self, member: member_type, id: int) -> None:
        super().__init__(member=member, id=id)
        self.reflection_matrix = get_r_of_fixed_end(m1=member)

    def get_equations(self) -> list:  # TODO add type
        # Assuming Wave equation to be A*e^(-I*k*x) + B*e^(-k*x) + C*e^(I*k*x) + D*e^(k*x)
        a1, b1, c1, d1 = self.member.get_parameters(id=self.id)

        eqns = self.reflection_matrix * Matrix([c1, d1]) - Matrix([a1, b1])

        return eqns
