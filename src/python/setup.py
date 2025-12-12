from pathlib import Path

from setuptools import find_packages, setup

ROOT = Path(__file__).resolve().parent

setup(
    name="ftst",
    version="0.1.0",
    description="Finite-temperature solver toolkit",
    packages=find_packages(),
    package_dir={"": "."},
    include_package_data=True,
    install_requires=[
        "numpy",
        "matplotlib",
        "pandas",
        "scipy",
    ],
)
