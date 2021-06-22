from setuptools import setup

setup(
    name = "dpamber",
    version = "0.0.5",
    install_requires = [
        "numpy",
        "dpdata>=0.2.1",
        "deepmd-kit",
        "tqdm",
    ],
    entry_points = {
        'console_scripts': ['dpamber=dpamber.cli:run'],
    },
    packages = ['dpamber'],
)
