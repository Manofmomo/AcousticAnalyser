from acoustic_analyser.modules.member import member as member_type
from acoustic_analyser.modules.rt_bc import get_r_of_fixed_end, get_r_of_free_end
from sympy import Matrix
import logging
from typing import List

logger = logging.getLogger("acoustic_analyser")


class bc:
    """This class emulates a boundary condition that can be added to members"""

    def __init__(self, member: member_type, id: int) -> None:
        if type(member) != member_type:
            raise Exception("members must belong to class Member")
        if not member.check_constraint_count():
            raise Exception(f"Member {member.id} already has 2 joints/BCs")
        member.increment_constraint_count()
        member.add_constraint(id=id)
        self.members: List[member_type] = [member]
        self.id = id
        self.theta = 0


class free_end(bc):
    """This is a free end boundary condition"""

    def __init__(self, member: member_type, id: int) -> None:
        super().__init__(member=member, id=id)
        self.reflection_matrix = get_r_of_free_end(m1=member)

        logger.debug(f"Reflection Matrix for free_end {self.id} calculated")

    def get_equations(self, w: float) -> list:
        a_plus, a_minus = self.members[0].get_parameters(id=self.id, w=w)

        matrix_reflect = self.reflection_matrix * a_plus - a_minus

        eqns = matrix_reflect
        logger.debug(f"Equations for free_end {self.id} calculated")
        return eqns


class fixed_end(bc):
    def __init__(self, member: member_type, id: int) -> None:
        super().__init__(member=member, id=id)
        self.reflection_matrix = get_r_of_fixed_end(m1=member)

        logger.debug(f"Reflection Matrix for fixed_end {self.id} calculated")

    def get_equations(self, w: float) -> list:
        """Gets the equations from the reflection and transmission matrices"""
        a_plus, a_minus = self.members[0].get_parameters(id=self.id, w=w)

        matrix_reflect = self.reflection_matrix * a_plus - a_minus

        eqns = matrix_reflect
        logger.debug(f"Equations for fixed_end {self.id} calculated")
        return eqns
