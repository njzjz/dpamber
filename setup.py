from setuptools import setup

setup(
    name = "dpamber",
    version = "0.0.8",
    install_requires = [
        "numpy",
        "dpdata[amber]>=0.2.1",
        "deepmd-kit",
        "tqdm",
    ],
    entry_points = {
        'console_scripts': ['dpamber=dpamber.cli:run'],
    },
    packages = ['dpamber'],
)
