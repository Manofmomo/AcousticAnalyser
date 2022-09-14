from src.acoustic_analyser import structure

test_structure = structure()

m1, _ = test_structure.add_member(length=0.5, rho=1, cross_section_area=1, youngs_modulus=1, inertia=1)
m2, _ = test_structure.add_member(length=0.5, rho=1, cross_section_area=1, youngs_modulus=1, inertia=1)
m3, _ = test_structure.add_member(length=1, rho=1, cross_section_area=1, youngs_modulus=1, inertia=1)

j, _ = test_structure.fixed_end(member=m1)
j1,_ =test_structure.two_member_joint(member_1=m1, member_2=m2, theta=0)
j2,_ =test_structure.two_member_joint(member_1=m2, member_2=m3, theta=0)

_ = test_structure.free_end(member=m3)

lhs=test_structure.get_equation_matrix()
print(lhs)