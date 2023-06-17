from pathlib import Path

import numpy as np

from dpamber.model_devi import calculate_devi


def test_model_devi(mocker):
    class FakeDP:
        def __init__(self, model) -> None:
            pass

        def eval(self, coords, cells, atom_types, mixed_type=False):
            return 0, np.zeros_like(coords), np.zeros(9)

    mocker.patch("deepmd.infer.DeepPot", FakeDP)

    # spy.assert_called_once_with(21)

    prefix = str(Path(__file__).parent / "corr/rc")
    calculate_devi(["graph.0.pb", "graph.1.pb"], 6, "corr/qmmm.parm7", ":1", prefix)
    x = np.loadtxt("model_devi.out", ndmin=2)
    assert x.shape == (1, 8)
