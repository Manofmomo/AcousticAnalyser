from sympy import symbols


class member:
    def __init__(
        self,
        length: float,
        rho: float,
        cross_section_area: float,
        youngs_modulus: float,
        inertia: float,
        omega: symbols,
        member_number: int,
    ) -> None:
        if length <= 0:
            raise ValueError("length must be greater than 0")
        if rho <= 0:
            raise ValueError("rho must be greater than 0")
        if cross_section_area <= 0:
            raise ValueError("cross_section_area must be greater than 0")
        if youngs_modulus <= 0:
            raise ValueError("youngs_modulus must be greater than 0")
        if inertia <= 0:
            raise ValueError("inertia must be greater than 0")

        self.length = length
        self.rho = rho
        self.cross_section_area = cross_section_area
        self.youngs_modulus = youngs_modulus
        self.inertia = inertia

        self.member_number = member_number
        self.omega = omega

        self.constraint_count = 0
        self.constraints = {}

        self.set_parameters()

    def check_constraint_count(self) -> bool:
        if self.constraint_count < 2:
            return True
        else:
            return False

    def increment_constraint_count(self):
        self.constraint_count = self.constraint_count + 1

    def add_constraint(self, id):
        self.constraints[id] = self.constraint_count

    def set_parameters(self) -> None:
        self.A, self.B, self.C, self.D = symbols(
            "A{i}, B{i}, C{i}, D{i}".format(i=self.member_number)
        )

    def get_parameters(self, id) -> list:
        if self.constraints[id]==1:
            return [self.A, self.B, self.C, self.D]
        else:
            return [self.C, self.D, self.A, self.B]