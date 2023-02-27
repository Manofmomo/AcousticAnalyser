from sympy import Matrix
from pickle import load as pickle_load
from acoustic_analyser.modules.member import member as member_type
import logging
import os

base_path = os.path.join(os.path.dirname(__file__), "..", "equations")

logger = logging.getLogger("acoustic_analyser")

file = open(os.path.join(base_path, "fixed_end.pkl"), mode="rb")
reflection_fixed = pickle_load(file)
file.close()

file = open(os.path.join(base_path, "free_end.pkl"), mode="rb")
reflection_free = pickle_load(file)
file.close()
logger.debug("Boundary Condition Equation Files Loaded")


def get_r_of_free_end(m1: member_type) -> Matrix:
    return reflection_free


def get_r_of_fixed_end(m1: member_type) -> Matrix:
    return reflection_fixed
