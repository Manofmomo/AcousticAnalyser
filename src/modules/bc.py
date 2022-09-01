from member import member as member_type


class bc:
    def __init__(self, member: member_type) -> None:
        if type(member) != member_type:
            raise Exception("members must belong to class Member")
        if not member.check_joint_count():
            raise Exception(f"Member {member.member_number} already has 2 joints/BCs")
        member.increment_joint_count()
        self.member = member


class free_end(bc):
    def __init__(self, member: member_type) -> None:
        super().__init__(member=member)

    def get_equations_from_members(self) -> tuple:
        self.v = self.members[0].get_equation()

    def get_equations(self) -> list:  # TODO add type
        self.get_equations_from_members()

        eqns = [None] * 2
        eqns[0] = self.strain_condition()
        eqns[1] = self.BM_condition()

        return eqns


class fixed_end(bc):
    def __init__(self, member: member_type) -> None:
        super().__init__(member=member)

    def get_equations_from_members(self) -> tuple:
        self.v = self.members[0].get_equation()

    def get_equations(self) -> list:  # TODO add type
        self.get_equations_from_members()

        eqns = [None] * 2
        eqns[0] = self.displacement_condition()
        eqns[1] = self.slope_condition()

        return eqns
