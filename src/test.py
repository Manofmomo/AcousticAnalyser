from acoustic_analyser import frame
import numpy as np
from scipy.optimize import newton

test_frame = frame.from_file(
    member_file="test_member.json", constraint_file="test_constraint.csv", debug=False
)

free_end = test_frame.fixed_end(member_id=0)
free_end = test_frame.free_end(member_id=1)
natural_freq = lambda x: test_frame.get_determinant(w=np.abs(x) * 2 * np.pi)
print(np.abs(newton(natural_freq, 0.1, tol=1e-09, maxiter=100)))
