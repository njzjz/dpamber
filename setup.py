from setuptools import setup

setup(
    name = "dpamber",
    version = "0.1.3",
    install_requires = [
        "numpy",
        "dpdata[amber]>=0.2.2",
        "deepmd-kit",
        "tqdm",
        "parmed",
        "ase",
        "scipy",
    ],
    entry_points = {
        'console_scripts': ['dpamber=dpamber.cli:run'],
        'dpdata.plugins': [
            'amber_md_qmmm=dpamber.read_amber:AmberMDQMMMFormat'
        ]
    },
    packages = ['dpamber'],
)
