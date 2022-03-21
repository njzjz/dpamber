import dpdata
import numpy as np
from dpdata.amber.mask import pick_by_amber_mask
from ase.geometry import wrap_positions, Cell
from typing import Union


def get_amber_fp(cutoff: float,
                 parmfile: str,
                 ncfile: str,
                 ll: str,
                 hl: str,
                 target: str = ":1",
                 out: str = None,
                 idx: Union[slice, list, int] = None,
                 suffix_mdfrc = None,
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

    Returns
    -------
    ms: dpdata.MultiSystems
        The output MultiSystems
    """
    ms = dpdata.MultiSystems()
    ep = r'@%EP'
    if cutoff > 0.:
        interactwith = "(%s)<@%f&!%s" % (target, cutoff, ep)
    else:
        interactwith = target

    if suffix_mdfrc is not None:
        ll_frc = ll + suffix_mdfrc
        hl_frc = hl + suffix_mdfrc
    else:
        ll_frc = None
        hl_frc = None

    s_ll = dpdata.LabeledSystem(
        ll, nc_file=ncfile, mdfrc_file=ll_frc, parm7_file=parmfile, fmt='amber/md/qmmm', qm_region=target)
    s_hl = dpdata.LabeledSystem(
        hl, nc_file=ncfile, mdfrc_file=hl_frc, parm7_file=parmfile, fmt='amber/md/qmmm', qm_region=target)
    if idx is not None:
        s_ll = s_ll[idx]
        s_hl = s_hl[idx]

    s_corr = s_ll.correction(s_hl)
    # wrap the coords...
    qm_index = pick_by_amber_mask(parmfile, target)
    cell = Cell(s_corr['cells'][0])
    wraped_coords = wrap_positions(cell.scaled_positions(s_corr['coords'][0]), cell=s_corr['cells'][0], pbc=True, center=cell.scaled_positions(np.mean(s_corr['coords'][0, qm_index], axis=0)))
    s_corr['coords'][0, :, :] = cell.cartesian_positions(wraped_coords)

    s_corr = s_corr.pick_by_amber_mask(
        parmfile, interactwith, pass_coords=True, nopbc=True)
    for ss in s_corr:
        ss = ss.remove_atom_names('EP')
        ms.append(ss)

    if out:
        ms.to_deepmd_npy(out)
    return ms


def run(args):
    get_amber_fp(cutoff=args.cutoff,
                 parmfile=args.parm7_file,
                 ncfile=args.nc,
                 ll=args.ll,
                 hl=args.hl,
                 target=args.qm_region,
                 out=args.out,
                 )
