from sympy import Matrix
from pickle import load as pickle_load
from src.modules.member import member as member_type
import logging

logger = logging.getLogger("acoustic_analyser")

file = open("src/equations/fixed_end", mode="rb")
reflection_fixed = pickle_load(file)
file.close()

file = open("src/equations/free_end", mode="rb")
reflection_free = pickle_load(file)
file.close()
logger.debug("Boundary Condition Equation Files Loaded")


def get_r_of_free_end(m1: member_type) -> Matrix:
    return reflection_free


def get_r_of_fixed_end(m1: member_type) -> Matrix:
    return reflection_fixed
