"""Calculate model deviation of the Amber MD traj."""
import numpy as np
import dpdata
from dpdata.amber.mask import load_param_file
from tqdm import tqdm


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
    import deepmd.DeepPot as DeepPot
    dps = [DeepPot(mm) for mm in models]

    sy = dpdata.System("rc", parm7_file=parm7_file,
                       fmt="amber/md", use_element_symbols=":1")
    parm7 = load_param_file(parm7_file)
    interactwith = "(%s)<@%f&!%s" % (qm_region, cutoff, r"@%EP")

    stds = []
    for ii, ssy in enumerate(tqdm(sy)):
        ssy = ssy.pick_by_amber_mask(
            parm7, maskstr=interactwith, pass_coords=True)
        ssy = ssy[0]
        ssy = ssy.remove_atom_names("EP")

        results = [ssy.predict(dp).pick_by_amber_mask(
            parm7, qm_region) for dp in dps]
        energies = [dp['energies'] for dp in results]
        energies = np.array(list(energies))
        forces = [dp['forces'] for dp in results]
        forces = np.array(list(forces))
        std_e_min = np.min(np.linalg.norm(np.std(energies, axis=0), axis=-1), axis=-1)
        std_e_ave = np.mean(np.linalg.norm(np.std(energies, axis=0), axis=-1), axis=-1)
        std_e_max = np.max(np.linalg.norm(np.std(energies, axis=0), axis=-1), axis=-1)
        std_f_min = np.min(np.linalg.norm(np.std(forces, axis=0), axis=-1), axis=-1)
        std_f_ave = np.mean(np.linalg.norm(np.std(forces, axis=0), axis=-1), axis=-1)
        std_f_max = np.max(np.linalg.norm(np.std(forces, axis=0), axis=-1), axis=-1)
        stds.append([ii, float(std_e_max), float(std_e_min), float(std_e_ave),  float(std_f_max), float(std_f_min), float(std_f_ave)])
    stds = np.array(stds)
    np.savetxt("model_devi.out", stds, header="step max_devi_e min_devi_e avg_devi_e max_devi_f min_devi_f avg_devi_f")

def run(args):
    calculate_devi(
        models = args.models,
        cutoff = args.cutoff,
        parm7_file = args.parm7_file,
        qm_region = args.qm_region
    )
