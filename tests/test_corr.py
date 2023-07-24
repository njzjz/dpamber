import shutil
from pathlib import Path

import dpdata
import numpy as np
import pytest

from dpamber.corr import get_amber_fp


@pytest.fixture
def test_system() -> dpdata.MultiSystems:
    """Tested system."""
    return dpdata.MultiSystems().from_deepmd_raw(
        Path(__file__).parent / "corr" / "dataset"
    )


def system_is_equal(sys1: dpdata.LabeledSystem, sys2: dpdata.LabeledSystem):
    """Check if two systems are equal.

    Parameters
    ----------
    sys1: dpdata.LabeledSystem
        The first system.
    sys2: dpdata.LabeledSystem
        The second system.

    Raises
    ------
    AssertionError
        If the two systems are not equal.
    """
    assert sys1["atom_names"] == sys2["atom_names"]
    assert sys1["atom_numbs"] == sys2["atom_numbs"]
    assert sys1.nopbc == sys2.nopbc
    np.testing.assert_almost_equal(sys1["atom_types"], sys2["atom_types"])
    np.testing.assert_almost_equal(sys1["coords"], sys2["coords"])
    np.testing.assert_almost_equal(sys1["cells"], sys2["cells"])
    np.testing.assert_almost_equal(sys1["energies"], sys2["energies"])
    np.testing.assert_almost_equal(sys1["forces"], sys2["forces"])
    np.testing.assert_almost_equal(sys1["aparam"], sys2["aparam"])


def get_single_system(multi: dpdata.MultiSystems) -> dpdata.LabeledSystem:
    """Assume the multi-system contains only one system.

    Parameters
    ----------
    multi: dpdata.MultiSystems
        The multi-system.

    Returns
    -------
    dpdata.LabeledSystem
        The system.
    """
    assert len(multi.systems) == 1
    return list(multi.systems.values())[0]


def test_corr(test_system):
    get_amber_fp(
        cutoff=6.0,
        parmfile=str(Path(__file__).parent / "corr/qmmm.parm7"),
        ncfile=str(Path(__file__).parent / "corr/rc.nc"),
        ll=str(Path(__file__).parent / "corr/high_level"),
        hl=str(Path(__file__).parent / "corr/low_level"),
        target=":1",
        out="tmp_data",
    )
    tmp_system = dpdata.MultiSystems().from_deepmd_npy("tmp_data")
    tmp_system_single = get_single_system(tmp_system)
    test_system_single = get_single_system(test_system)
    system_is_equal(tmp_system_single, test_system_single)


def test_corr_hdf5(test_system):
    get_amber_fp(
        cutoff=6.0,
        parmfile=str(Path(__file__).parent / "corr/qmmm.parm7"),
        ncfile=str(Path(__file__).parent / "corr/rc.nc"),
        ll=str(Path(__file__).parent / "corr/high_level"),
        hl=str(Path(__file__).parent / "corr/low_level"),
        target=":1",
        out="tmp_data.hdf5",
    )
    tmp_system = dpdata.MultiSystems().from_deepmd_hdf5("tmp_data.hdf5")
    tmp_system_single = get_single_system(tmp_system)
    test_system_single = get_single_system(test_system)
    system_is_equal(tmp_system_single, test_system_single)


def test_md_corr(test_system):
    shutil.copyfile(
        Path(__file__).parent / "corr/low_level.mdfrc",
        Path(__file__).parent / "md/md.mdfrc",
    )
    get_amber_fp(
        cutoff=6.0,
        parmfile=str(Path(__file__).parent / "corr/qmmm.parm7"),
        ncfile=str(Path(__file__).parent / "corr/rc.nc"),
        ll=str(Path(__file__).parent / "corr/high_level"),
        hl=str(Path(__file__).parent / "md/md"),
        target=":1",
        out="tmp_data",
    )
    tmp_system = dpdata.MultiSystems().from_deepmd_npy("tmp_data")
    tmp_system_single = get_single_system(tmp_system)
    test_system_single = get_single_system(test_system)
    system_is_equal(tmp_system_single, test_system_single)


def test_uncoverage_system():
    tmp_system = dpdata.LabeledSystem(
        str(Path(__file__).parent / "uncoverage/uncoverage"),
        nc_file=str(Path(__file__).parent / "corr/rc.nc"),
        mdfrc_file=str(Path(__file__).parent / "corr/low_level.mdfrc"),
        parm7_file=str(Path(__file__).parent / "corr/qmmm.parm7"),
        fmt="amber/md/qmmm",
        qm_region=":1",
    )
    assert len(tmp_system) == 0


def test_uncoverage_corr():
    shutil.copyfile(
        Path(__file__).parent / "corr/low_level.mdfrc",
        Path(__file__).parent / "uncoverage/uncoverage.mdfrc",
    )
    tmp_system = get_amber_fp(
        cutoff=6.0,
        parmfile=str(Path(__file__).parent / "corr/qmmm.parm7"),
        ncfile=str(Path(__file__).parent / "corr/rc.nc"),
        ll=str(Path(__file__).parent / "corr/high_level"),
        hl=str(Path(__file__).parent / "uncoverage/uncoverage"),
        target=":1",
    )
    assert len(tmp_system) == 0
