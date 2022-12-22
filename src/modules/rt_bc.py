from sympy import Matrix
from pickle import load as pickle_load
from src.modules.member import member as member_type

file = open("src/equations/fixed_end", mode="rb")
reflection_fixed = pickle_load(file)
file.close()

file = open("src/equations/free_end", mode="rb")
reflection_free = pickle_load(file)
file.close()

def get_r_of_free_end(m1: member_type) -> Matrix:
    return reflection_free


def get_r_of_fixed_end(m1: member_type) -> Matrix:
    return reflection_fixed
