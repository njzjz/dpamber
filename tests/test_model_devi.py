from pathlib import Path

import numpy as np

from dpamber.model_devi import calculate_devi


def test_model_devi(mocker):
    class FakeDP:
        def __init__(self, model) -> None:
            pass

        def eval(self, coords, cells, atom_types, mixed_type=False):
            return 0, np.zeros_like(coords), np.zeros(9)

        def get_type_map(self):
            return ["C", "H", "O", "N", "P", "S", "HW", "OW"]

    mocker.patch("deepmd.infer.DeepPot", FakeDP)

    # spy.assert_called_once_with(21)

    calculate_devi(
        ["graph.0.pb", "graph.1.pb"],
        6,
        str(Path(__file__).parent / "corr/qmmm.parm7"),
        ":1",
        str(Path(__file__).parent / "corr/rc"),
    )
    x = np.loadtxt("model_devi.out", ndmin=2)
    assert x.shape == (1, 8)
