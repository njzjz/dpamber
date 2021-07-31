import dpdata
import numpy as np


def get_amber_fp(cutoff: float,
                 parmfile: str,
                 ncfile: str,
                 ll: str,
                 hl: str,
                 target: str = ":1",
                 out: str = None):
    ms = dpdata.MultiSystems()
    ep = r'@%EP'
    if cutoff > 0.:
        interactwith = "(%s)<@%f&!%s" % (target, cutoff, ep)
    else:
        interactwith = target

    s_ll = dpdata.LabeledSystem(
        ll, nc_file=ncfile, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
    s_hl = dpdata.LabeledSystem(
        hl, nc_file=ncfile, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
    s_corr = s_ll.correction(s_hl)
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
