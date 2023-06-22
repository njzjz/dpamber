import numpy as np
from dpdata.system import Axis, DataType, LabeledSystem, System

System.DTYPES = System.DTYPES + (
    DataType("aparam", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 1), required=False),
)
LabeledSystem.DTYPES = LabeledSystem.DTYPES + (
    DataType("aparam", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 1), required=False),
)
ep = None
