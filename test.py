#%%
from src.acoustic_analyser import frame

test_frame = frame.from_file(
    member_file="test_member.json", constraint_file="test_constraint.csv"
)

free_end = test_frame.fixed_end(member_id=0)
free_end = test_frame.free_end(member_id=1)

#%%
eqs = test_frame.get_equation_matrix()


#%%
test_frame.constraints[0].reflection_matrix_11.expand()