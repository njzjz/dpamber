import dpdata
import numpy as np


def get_amber_fp(cutoff: float,
                 target: str = ":1",
                 out: str = None,
            ) -> dpdata.MultiSystems:
    """(qmwater-qmwater_not)-(mmwater-mmwater_not)
    = (qmwater-mmwater) - (qmwater_not-mmwater_not)
    """
    ms = dpdata.MultiSystems()
    ep = r'@%EP'
    if cutoff > 0.:
        interactwith = "(%s)<@%f&!%s" % (target, cutoff, ep)
    else:
        interactwith = target

    s_qmwater = dpdata.LabeledSystem(
        "qmwater", fmt='amber/md', use_element_symbols=target)
    s_qmwater = s_qmwater.remove_atom_names('EP')
    s_mmwater = dpdata.LabeledSystem(
        "mmwater", fmt='amber/md', use_element_symbols=target)
    s_mmwater = s_mmwater.remove_atom_names('EP')
    assert(np.all(s_qmwater['coords']==s_mmwater['coords']))
    s_allcorr = s_mmwater.correction(s_qmwater)
    
    s_qmwaternot = dpdata.LabeledSystem(
        "qmwater_not", fmt='amber/md')
    s_qmwaternot = s_qmwaternot.remove_atom_names('EP')
    s_mmwaternot = dpdata.LabeledSystem(
        "mmwater_not", fmt='amber/md')
    s_mmwaternot = s_mmwaternot.remove_atom_names('EP')
    assert(np.all(s_qmwaternot['coords']==s_mmwaternot['coords']))
    s_notcorr = s_mmwaternot.correction(s_qmwaternot)

    # manually process
    natom_not = s_notcorr.get_natoms()
    assert(np.all(s_allcorr['coords'][:,-natom_not:]==s_notcorr['coords']))
    s_allcorr.data['energies'] -= s_notcorr['energies']
    s_allcorr.data['forces'][:,-natom_not:] -= s_notcorr['forces']

    s_allcorr = s_allcorr.pick_by_amber_mask(
        "noepw.parm7", interactwith, pass_coords=True, nopbc=True)
    for ss in s_allcorr:
        ms.append(ss)

    if out:
        ms.to_deepmd_npy(out)
    return ms


def run(args):
    get_amber_fp(cutoff=args.cutoff,
                 target=args.qm_region,
                 out=args.out,
                 )
