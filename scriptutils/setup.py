from setuptools import find_packages, setup

setup(
    name="scriptutils",
    version=0.1,
    author="Lukas Nickel",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "astropy<6",
        "numpy",
        "pydantic<2.0",
        "rich",
    ],
    license="MIT",
)
