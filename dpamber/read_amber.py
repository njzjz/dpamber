import re
import os
from scipy.io import netcdf
import numpy as np
from dpdata.amber.mask import pick_by_amber_mask
from dpdata.unit import EnergyConversion
from dpdata.periodic_table import ELEMENTS
from dpdata.format import Format
from dpdata.plugins.amber import AmberMDFormat
from ase.geometry import cellpar_to_cell


kcalmol2eV = EnergyConversion("kcal_mol", "eV").value()
symbols = ['X'] + ELEMENTS

energy_convert = kcalmol2eV
force_convert = energy_convert


def read_amber_traj(parm7_file, nc_file, mdfrc_file=None, mden_file=None, mdout_file=None,
                    qm_region=None, labeled=True,
                    ):
    """The amber trajectory includes:
    * nc, NetCDF format, stores coordinates
    * mdfrc, NetCDF format, stores forces
    * mden (optional), text format, stores energies
    * mdout (optional), text format, may store energies if there is no mden_file
    * parm7, text format, stores types

    Parameters
    ----------
    parm7_file, nc_file, mdfrc_file, mden_file, mdout_file: str
      filenames
    qm_region: None or list or str
      amber mask
    """

    flag_atom_type = False
    flag_atom_numb = False
    amber_types = []
    atomic_number = []
    with open(parm7_file) as f:
        for line in f:
            if line.startswith("%FLAG"):
                flag_atom_type = line.startswith("%FLAG AMBER_ATOM_TYPE")
                flag_atom_numb = line.startswith("%FLAG ATOMIC_NUMBER")
            elif flag_atom_type or flag_atom_numb:
                if line.startswith("%FORMAT"):
                    fmt = re.findall(r'\d+', line)
                    fmt0 = int(fmt[0])
                    fmt1 = int(fmt[1])
                else:
                    for ii in range(fmt0):
                        start_index = ii * fmt1
                        end_index = (ii + 1) * fmt1
                        if end_index >= len(line):
                            continue
                        content = line[start_index:end_index].strip()
                        if flag_atom_type:
                            amber_types.append(content)
                        elif flag_atom_numb:
                            atomic_number.append(int(content))
    ml_types = []
    if qm_region is None:
        qm_region = []

    if isinstance(qm_region, str):
        qm_region = pick_by_amber_mask(
            parm7_file, qm_region)
    
    for ii, (tt, nn) in enumerate(zip(amber_types, atomic_number)):
        if ii in qm_region:
            ml_types.append(symbols[nn])
        elif tt in ("OW", "HW", "EP"):
            ml_types.append(tt)
        else:
            ml_types.append("m" + symbols[nn])

    with netcdf.netcdf_file(nc_file, 'r') as f:
        coords = np.array(f.variables["coordinates"][:])
        cell_lengths = np.array(f.variables["cell_lengths"][:])
        cell_angles = np.array(f.variables["cell_angles"][:])
        shape = cell_lengths.shape
        cells = np.zeros((shape[0], 3, 3))
        for ii in range(cell_lengths.shape[0]):
            cells[ii, :, :] = cellpar_to_cell([*cell_lengths[ii], *cell_angles[ii]])

    if labeled:
        with netcdf.netcdf_file(mdfrc_file, 'r') as f:
            forces = np.array(f.variables["forces"][:])

        # load energy from mden_file or mdout_file
        energies = []
        if mden_file is not None and os.path.isfile(mden_file):
            with open(mden_file) as f:
                for line in f:
                    if line.startswith("L6"):
                        s = line.split()
                        if s[2] != "E_pot":
                            energies.append(float(s[2]))
        elif mdout_file is not None and os.path.isfile(mdout_file):
            with open(mdout_file) as f:
                for line in f:
                    if "EPtot" in line:
                        s = line.split()
                        energies.append(float(s[-1]))
        else:
            raise RuntimeError(
                "Please provide one of mden_file and mdout_file")

    atom_names, atom_types, atom_numbs = np.unique(
        ml_types, return_inverse=True, return_counts=True)

    data = {}
    data['atom_names'] = list(atom_names)
    data['atom_numbs'] = list(atom_numbs)
    data['atom_types'] = atom_types
    if labeled:
        data['forces'] = forces * force_convert
        data['energies'] = np.array(energies) * energy_convert
    data['coords'] = coords
    data['cells'] = cells
    data['orig'] = np.array([0, 0, 0])
    return data


@Format.register("amber/md/qmmm")
class AmberMDQMMMFormat(AmberMDFormat):
    def from_system(self, file_name=None, parm7_file=None, nc_file=None, qm_region=None, **kwargs):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        return read_amber_traj(parm7_file=parm7_file, nc_file=nc_file, qm_region=qm_region, labeled=False)

    def from_labeled_system(self, file_name=None, parm7_file=None, nc_file=None, mdfrc_file=None, mden_file=None, mdout_file=None, qm_region=None, **kwargs):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        if mdfrc_file is None:
            mdfrc_file = file_name + ".mdfrc"
        if mden_file is None:
            mden_file = file_name + ".mden"
        if mdout_file is None:
            mdout_file = file_name + ".mdout"
        return read_amber_traj(parm7_file, nc_file, mdfrc_file, mden_file, mdout_file, qm_region)
