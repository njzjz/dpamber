import dpdata
import numpy as np


def get_amber_fp(cutoff: float,
                 parmfile: str,
                 ll: str,
                 hl: str,
                 target: str = ":1",
                 out: str = "dataset"):
    ms = dpdata.MultiSystems()
    ep = r'@%EP'
    interactwith = "(%s)<:%f&!%s" % (target, cutoff, ep)

    s_ll = dpdata.LabeledSystem(
        ll, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
    s_hl = dpdata.LabeledSystem(
        hl, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
    s_corr = s_ll.correction(s_hl)
    s_corr = s_corr.pick_by_amber_mask(
        parmfile, interactwith, pass_coords=True)
    for ss in s_corr:
        ss = ss.remove_atom_names('EP')
        ms.append(ss)

    ms.to_deepmd_npy(out)


def run(args):
    get_amber_fp(cutoff=args.cutoff,
                 parmfile=args.parm7_file,
                 ll=args.ll,
                 hl=args.hl,
                 target=args.qm_region,
                 out=args.out,
                 )
