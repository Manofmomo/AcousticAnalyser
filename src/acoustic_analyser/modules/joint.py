from asyncio.log import logger
from acoustic_analyser.modules.member import member as member_type
from acoustic_analyser.modules.rt_joint import get_rt_of_two_member
from sympy import Matrix
from typing import List
from math import pi
import logging

logger = logging.getLogger("acoustic_analyser")


class joint:
    """This class emulates a joint that can be added to members"""

    def __init__(self, members: List[member_type], id: int) -> None:
        for member in members:
            if type(member) != member_type:
                raise Exception("members must belong to class Member")
            if not member.check_constraint_count():
                raise Exception(f"Member {member.id} already has 2 joints/BCs")
        for member in members:
            member.increment_constraint_count()
            member.add_constraint(id=id)

        self.members: List[member_type] = members
        self.id = id


class two_member(joint):
    """This class emulates a two-member joint"""

    def __init__(
        self, theta: float, member_1: member_type, member_2: member_type, id: int
    ) -> None:
        self.theta = theta * pi / 180  # Degree to radians
        self.member_joint_type = []
        members = [member_1, member_2]
        super().__init__(members=members, id=id)

    def set_rt_matrices(self, w: float):
        (
            self.reflection_matrix_11,
            self.reflection_matrix_22,
            self.transmission_matrix_12,
            self.transmission_matrix_21,
        ) = get_rt_of_two_member(
            m1=self.members[0], m2=self.members[1], theta=self.theta, w=w
        )

    def get_equations(self, w: float) -> list:
        """Gets the equations from the reflection and transmission matrices"""
        a_plus, a_minus = self.members[0].get_parameters(id=self.id, w=w)
        b_minus, b_plus = self.members[1].get_parameters(id=self.id, w=w)

        self.set_rt_matrices(w=w)

        matrix_reflect = (
            self.reflection_matrix_11 * a_plus
            + self.transmission_matrix_21 * b_minus
            - a_minus
        )
        # DIFFERENT REFLECTION TRANSMISSION
        matrix_transmit = (
            self.transmission_matrix_12 * a_plus
            + self.reflection_matrix_22 * b_minus
            - b_plus
        )
        eqns = matrix_reflect.col_join(matrix_transmit)
        logger.debug(f"Equations for joint {self.id} calculated")
        return eqns
