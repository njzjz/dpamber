from setuptools import setup

setup(
    name = "dpamber",
    version = "0.0.1",
    install_requires = [
        "numpy",
        "dpdata",
        "deepmd-kit",
    ],
    entry_points = {
        'console_scripts': ['dpamber=dpamber.cli:run'],
    },
    packages = ['dpamber'],
)