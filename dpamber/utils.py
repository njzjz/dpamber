import dpdata
import numpy as np
from ase.geometry import Cell, get_distances, wrap_positions
from dpdata.amber.mask import pick_by_amber_mask


def wrap_system(
    s_corr: dpdata.System, parmfile: str, target: str, cutoff: float
) -> dpdata.MultiSystems:
    """Wrap the coordinates around the QM region.

    Parameters
    ----------
    s_corr: dpdata.System
        The system to be processed.
    parmfile: str
        The parm7 file.
    target: str
        The QM system mask.
    cutoff: float
        The QM/MM cutoff radius.

    Returns
    -------
    ms: dpdata.MultiSystems
        The output MultiSystems
    """
    ms = dpdata.MultiSystems()
    ep = r"@%EP"
    if cutoff > 0.0:
        interactwith = f"({target})<@{cutoff:f}&!{ep}"
    else:
        interactwith = target
    # wrap the coords...
    qm_index = pick_by_amber_mask(parmfile, target)
    qm_coords = s_corr["coords"][:, qm_index, :]
    for ii in range(len(s_corr)):
        cell = Cell(s_corr["cells"][ii])
        qm_coord = qm_coords[ii]
        qm_distances = get_distances(qm_coord, cell=cell, pbc=True)[1]
        # find the coord that has the minimal total distance as the center
        center = qm_coords[ii, np.argmin(np.sum(qm_distances, axis=1))]
        wraped_coords = wrap_positions(
            s_corr["coords"][ii],
            cell=s_corr["cells"][ii],
            pbc=True,
            center=cell.scaled_positions(center),
        )
        s_corr["coords"][ii, :, :] = wraped_coords

    s_corr = s_corr.pick_by_amber_mask(
        parmfile, interactwith, pass_coords=True, nopbc=True
    )
    for ss in s_corr:
        if "EP" in ss["atom_names"]:
            ss = ss.remove_atom_names("EP")
        ms.append(ss)
    return ms
