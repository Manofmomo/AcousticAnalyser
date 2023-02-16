from acoustic_analyser.acoustic_analyser import frame
import numpy as np
from sympy import Matrix

test_frame = frame()

test_frame.add_member(
    length=0.5, density=7800, youngs_modulus=206e9, height=1.27e-2, id=0
)
test_frame.add_member(
    length=0.5, density=7800, youngs_modulus=206e9, height=1.27e-2, id=1
)

free_end = test_frame.fixed_end(member_id=0)
joint = test_frame.two_member_joint(member_1_id=0, member_2_id=1, theta=45)
free_end = test_frame.free_end(member_id=1)


matrix_1 = test_frame.get_equation_matrix(w=10)

print(abs(np.linalg.det(matrix_1)))


c1_plus, c1_minus = test_frame.members[0].get_parameters(w=3, id=1)
c2_minus, c2_plus = test_frame.members[1].get_parameters(w=3, id=1)
a_minus, a_plus = test_frame.members[0].get_parameters(w=3, id=0)
b_plus, b_minus = test_frame.members[1].get_parameters(w=3, id=2)
A = Matrix([])
# Fixed End Boundary Condition
A = A.col_join(test_frame.constraints[0].reflection_matrix * a_minus - a_plus)
# Joint
A = A.col_join(
    test_frame.constraints[1].reflection_matrix_11 * c1_plus
    + test_frame.constraints[1].transmission_matrix_21 * c2_minus
    - c1_minus
)
A = A.col_join(
    test_frame.constraints[1].reflection_matrix_22 * c2_minus
    + test_frame.constraints[1].transmission_matrix_12 * c1_plus
    - c2_plus
)
# Free End
A = A.col_join(test_frame.constraints[2].reflection_matrix * b_plus - b_minus)
# Propagation
A = A.col_join(test_frame.members[0].propagation_matrix_subs * a_plus - c1_plus)
A = A.col_join(test_frame.members[0].propagation_matrix_subs * c1_minus - a_minus)
A = A.col_join(test_frame.members[1].propagation_matrix_subs * c2_plus - b_plus)
A = A.col_join(test_frame.members[1].propagation_matrix_subs * b_minus - c2_minus)
A.expand() - test_frame.eqn_matrix.expand()
