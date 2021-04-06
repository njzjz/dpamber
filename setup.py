from setuptools import setup

setup(
    name = "dpamber",
    version = "0.0.3",
    install_requires = [
        "numpy",
        "dpdata",
        "deepmd-kit",
        "tqdm",
    ],
    entry_points = {
        'console_scripts': ['dpamber=dpamber.cli:run'],
    },
    packages = ['dpamber'],
)