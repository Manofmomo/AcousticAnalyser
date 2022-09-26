from sympy import symbols, Matrix


class member:
    # This class emulates a member with all its physical properties
    def __init__(
        self,
        length: float,
        rho: float,
        cross_section_area: float,
        youngs_modulus: float,
        inertia: float,
        omega: symbols,
        id: int,
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

        self.id = id
        self.omega = omega

        self.constraint_count = 0
        self.constraint_ids = []

        self.set_parameters()

    def check_constraint_count(self) -> bool:
        # Each member can have only 2 constraints added to it at present
        if self.constraint_count < 2:
            return True
        else:
            return False

    def increment_constraint_count(self):
        self.constraint_count = self.constraint_count + 1

    def add_constraint(self, id):
        self.constraint_ids.append(id)

    def set_parameters(self) -> None:
        a_b_plus, a_e_plus, a_b_minus, a_e_minus, a_l_plus, a_l_minus = symbols(
            "a_b^+{i}, a_e^+{i}, a_b^-{i}, a_e^-{i}, a_l^+{i}, a_l^-{i}".format(
                i=self.id
            )
        )
        self.a_plus = Matrix([a_b_plus, a_e_plus, a_l_plus])
        self.a_minus = Matrix([a_b_minus, a_e_minus, a_l_minus])

        self.b_plus = self.a_plus # Have to add propagation matrices
        self.b_minus = self.a_minus

    def get_parameters(self, id: int = None) -> list:
        # This function gives back the set of parameters to be used.
        # It corrects for the sign convention of the constraint when returning parameters
        # Positive direction is from lower to higher constraint id
        if id == max(self.constraint_ids):
            return [self.b_minus, self.b_plus]
        else:
            return [self.a_plus, self.a_minus]
