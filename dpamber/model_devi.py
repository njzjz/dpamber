"""Calculate model deviation of the Amber MD traj."""
import dpdata
import numpy as np
from dpdata.amber.mask import load_param_file
from tqdm import tqdm


def calculate_devi(
    models: list,
    cutoff: float,
    parm7_file: str = "qmmm.parm7",
    qm_region: str = ":1",
    prefix: str = "rc",
):
    """Calculate model deviation of the Amber MD traj.

    Parameters
    ----------
    models: list[str]
        models, for example ["graph.0.pb", "graph.1.pb"]
    cutoff: float
        The QM/MM cutoff radius.
    parm7_file: str
        The parm7 file.
    qm_region: str
        The QM system mask.
    prefix: str
        The prefix of MD files.
    """
    # Deep potentials
    from deepmd.infer import DeepPot, calc_model_devi

    dps = [DeepPot(mm) for mm in models]

    sy = dpdata.System(
        file_name=prefix,
        parm7_file=parm7_file,
        fmt="amber/md/qmmm",
        qm_region=qm_region,
    )
    parm7 = load_param_file(parm7_file)
    ep = r"@%EP"
    if cutoff > 0.0:
        interactwith = f"({qm_region})<@{cutoff:f}&!{ep}"
    else:
        interactwith = qm_region

    stds = []
    for ii, ss in enumerate(tqdm(sy)):
        ms = ss.pick_by_amber_mask(parm7, interactwith, pass_coords=True)
        ss = list(ms.systems.values())[0]
        if "EP" in ss["atom_names"]:
            ss = ss.remove_atom_names("EP")
        ss.apply_type_map(dps[0].get_type_map())
        devi = calc_model_devi(ss["coords"], ss["cells"], ss["atom_types"], dps)
        devi[0, 0] = ii
        stds.append(devi)
    stds = np.concatenate(stds)
    np.savetxt(
        "model_devi.out",
        stds,
        header="step max_devi_v min_devi_v avg_devi_v max_devi_f min_devi_f avg_devi_f devi_e",
    )


def run(args):
    calculate_devi(
        models=args.models,
        cutoff=args.cutoff,
        parm7_file=args.parm7_file,
        qm_region=args.qm_region,
        prefix=args.prefix,
    )
