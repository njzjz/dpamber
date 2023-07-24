import numpy as np
from dpdata.data_type import Axis, DataType, register_data_type

register_data_type(
    DataType("drdq", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 3, -1), required=False)
)
ep = None
