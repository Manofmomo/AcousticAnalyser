from src.acoustic_analyser import frame
import numpy as np

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


matrix_1 = test_frame.get_equation_matrix(w=3)
matrix_2 = test_frame.get_equation_matrix(w=10)

print(abs(np.linalg.det(matrix_1)))
print(abs(np.linalg.det(matrix_10)))
