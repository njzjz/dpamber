import dpdata
import numpy as np
#from multiprocessing import Pool

def read_traj(prefix, nc=None, frc_suffix=None):
    s_ll = dpdata.LabeledSystem(
        ll, mdfrc_file=ll+".frc.nc", fmt='amber/md', use_element_symbols="(%s)&!%s" % (target, ep))
    # remove EP in QM region
    s_ll = s_ll.pick_by_amber_mask(
        ll+".parm7", ":* & !((%s)&%s)" % (target, ep))


def get_amber_fp(cutoff: float,
                 parmfile: str,
                 ll: str,
                 hl: str,
                 type_map: list,
                 target: str = ":1",
                 out: str = "dataset",
                 nc: str = None,
                 frc_suffix: str = None):
    ms = dpdata.MultiSystems()
    ep = r'@%EP'
    interactwith = "(%s)<:%f&!%s" % (target, cutoff, ep)

    s_ll = dpdata.LabeledSystem(
        ll, mdfrc_file=ll+".frc.nc", fmt='amber/md', use_element_symbols="(%s)&!%s" % (target, ep))
    # remove EP in QM region
    s_ll = s_ll.pick_by_amber_mask(
        ll+".parm7", ":* & !((%s)&%s)" % (target, ep))
    s_hl = dpdata.LabeledSystem(
        hl, mdfrc_file=hl+".frc.nc", fmt='amber/md', use_element_symbols=target)
    s_corr = s_ll.correction(s_hl)
    s_corr = s_corr.pick_by_amber_mask(
        hl+".parm7", interactwith, pass_coords=True)
    for ss in s_corr:
        ss = ss.remove_atom_names('EP')
        ms.append(ss)

    ms.to_deepmd_npy(out, set_size=999999)

