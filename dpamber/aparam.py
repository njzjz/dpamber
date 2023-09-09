import numpy as np
from dpdata.data_type import Axis, DataType, register_data_type

register_data_type(
    DataType("aparam", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 1), required=False),
    labeled=False,
)
ep = None
