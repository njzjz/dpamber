from typing import List, Optional, Union

import dpdata
import numpy as np
from ase.geometry import Cell, get_distances, wrap_positions
from dpdata.amber.mask import pick_by_amber_mask

from dpamber.disang import Disang


def get_amber_fp(
    cutoff: float,
    parmfile: str,
    ncfile: str,
    ll: str,
    hl: str,
    target: str = ":1",
    out: Optional[str] = None,
    idx: Optional[Union[slice, list, int]] = None,
    suffix_mdfrc=None,
    disang_file: Optional[str] = None,
    rxn_idx: Optional[List[int]] = None,
) -> dpdata.MultiSystems:
    """Use Ambertools to do correction calculation between a high level potential and a low level potential.

    Parameters
    ----------
    cutoff: float
        The QM/MM cutoff radius.
    parmfile: str
        The original parm file.
    ncfile: str
        The coordinates file.
    ll: str
        The low level system prefix.
    hl: str
        The high level system prefix.
    target: str
        The QM system mask.
    out: str
        The output deepmd/npy directory.
    idx: slice or list or int
        index
    suffix_mdfrc: str, optional
        suffix of mdfrc file
    disang_file: str, optional
        The AMBER disang file for generalized forces
    rxn_idx: list of int, optional
        index of reaction coordinates

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

    if suffix_mdfrc is not None:
        ll_frc = ll + suffix_mdfrc
        hl_frc = hl + suffix_mdfrc
    else:
        ll_frc = None
        hl_frc = None

    if disang_file is not None:
        disang = Disang(disang_file)
        if rxn_idx is not None:
            disang.restraints = [
                rr for ii, rr in enumerate(disang.restraints) if ii in rxn_idx
            ]
    else:
        disang = None

    s_ll = dpdata.LabeledSystem(
        ll,
        nc_file=ncfile,
        mdfrc_file=ll_frc,
        parm7_file=parmfile,
        fmt="amber/md/qmmm",
        qm_region=target,
        exclude_unconverged=False,
        disang=disang,
    )
    s_hl = dpdata.LabeledSystem(
        hl,
        nc_file=ncfile,
        mdfrc_file=hl_frc,
        parm7_file=parmfile,
        fmt="amber/md/qmmm",
        qm_region=target,
        exclude_unconverged=False,
        disang=disang,
    )
    if idx is not None:
        s_ll = s_ll[idx]
        s_hl = s_hl[idx]

    s_corr = s_ll.correction(s_hl)

    # remove unconverged frames
    idx_pick = ~np.logical_or(
        np.isnan(s_corr["energies"]), np.isnan(s_corr["forces"]).any(axis=(1, 2))
    )
    s_corr = s_corr[idx_pick]

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
        for name, numb in zip(ss["atom_names"], ss["atom_numbs"]:
            if numb == 0:
                ss = ss.remove_atom_names(name)
        ms.append(ss)

    if out:
        if out.endswith(".hdf5"):
            ms.to_deepmd_hdf5(out)
        else:
            ms.to_deepmd_npy(out)
    return ms


def run(args):
    get_amber_fp(
        cutoff=args.cutoff,
        parmfile=args.parm7_file,
        ncfile=args.nc,
        ll=args.ll,
        hl=args.hl,
        target=args.qm_region,
        out=args.out,
        suffix_mdfrc=args.suffix_mdfrc,
        disang_file=args.disang,
        rxn_idx=args.rxn,
    )
