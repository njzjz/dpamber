"""Calculate model deviation of the Amber MD traj."""
import numpy as np
import deepmd.DeepPot as DeepPot
import dpdata
from dpdata.amber.mask import load_param_file


def calculate_devi(models: list,
                   cutoff: float,
                   parm7_file: str = "qmmm.parm7",
                   qm_region: str = ":1"
                   ):
    """
    Parameters
    ----------
    models: ["graph.0.pb", "graph.1.pb"]
    """
    # Deep potentials
    dps = [DeepPot(mm) for mm in models]

    sy = dpdata.System("rc", parm7_file=parm7_file,
                       fmt="amber/md", use_element_symbols=":1")
    parm7 = load_param_file(parm7_file)
    interactwith = "(%s)<:%f&!%s" % (qm_region, cutoff, r"@%EP")

    stds = []
    for ssy in sy:
        ssy = ssy.pick_by_amber_mask(
            parm7, maskstr=interactwith, pass_coords=True)
        ssy = ssy[0]
        ssy = ssy.remove_atom_names("EP")

        forces = [ssy.predict(dp).pick_by_amber_mask(
            parm7, qm_region)['forces'] for dp in dps]
        forces = np.array(list(forces))
        std = np.max(np.linalg.norm(np.std(forces, axis=0), axis=-1), axis=-1)
        stds.append(float(std))
    stds = np.array(stds)
    np.savetxt("std.txt", stds)

def run(args):
    calculate_devi(
        models = args.models,
        cutoff = args.cutoff,
        parm7_file = args.parm7_file,
        qm_region = args.qm_region
    )
