from src.modules.member import member as member_type
from src.modules.rt_bc import get_r_of_fixed_end, get_r_of_free_end
from sympy import Matrix
import logging

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s %(name)s %(levelname)s:%(message)s"
)
logger = logging.getLogger(__name__)


class bc:
    # This class emulates a boundary condition that can be added to members
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
    # This is a free end boundary condition
    def __init__(self, member: member_type, id: int) -> None:
        super().__init__(member=member, id=id)
        self.reflection_matrix = get_r_of_free_end(m1=member)
        print(self.reflection_matrix)
        logger.debug(f"Reflection Matrix for free_end {self.id} calculated")

    def get_equations(self) -> list:
        a_plus, a_minus = self.member.get_parameters(id=self.id)

        matrix_reflect = self.reflection_matrix * a_plus - a_minus

        eqns = matrix_reflect
        logger.debug(f"Equations for free_end {self.id} calculated")
        return eqns


class fixed_end(bc):
    def __init__(self, member: member_type, id: int) -> None:
        super().__init__(member=member, id=id)
        self.reflection_matrix = get_r_of_fixed_end(m1=member)
        print(self.reflection_matrix)
        logger.debug(f"Reflection Matrix for free_end {self.id} calculated")

    def get_equations(self) -> list:
        # Gets the equations from the reflection and transmission matrices
        a_plus, a_minus = self.member.get_parameters(id=self.id)

        matrix_reflect = self.reflection_matrix * a_plus - a_minus

        eqns = matrix_reflect
        logger.debug(f"Equations for fixed_end {self.id} calculated")
        return eqns
