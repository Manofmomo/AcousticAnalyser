from sympy import Matrix
from pickle import load as pickle_load
from src.modules.member import member as member_type


def get_r_of_free_end(m1: member_type) -> Matrix:
    file = open("src/equations/free_end", mode="rb")
    reflection = pickle_load(file)
    return reflection


def get_r_of_fixed_end(m1: member_type) -> Matrix:
    file = open("src/equations/fixed_end", mode="rb")
    reflection = pickle_load(file)
    return reflection
