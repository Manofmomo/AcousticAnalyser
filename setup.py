from setuptools import setup

setup(
    name="acoustic_analyser",
    version="0.0.9",
    author="Dhruv Bhatia, Dr. Abhijit Sarkar",
    description="An Open Source package to perform nodal analysis on simple 2D structures using wave based approach",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Manofmomo/AcousticAnalyser",
    package_dir={"": "src"},
    packages=["acoustic_analyser"],
    install_requires=[
        "sympy",
        "numpy",
        "scipy",
        "matplotlib",
    ],
    include_package_data=True,
)
