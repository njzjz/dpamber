import dpdata
from dpdata.amber.mask import load_param_file, pick_by_amber_mask
import parmed.tools as PT
from copy import copy


def qmbuffer(cutoff: float,
              ncfile: str,
              parmfile: str,
              hl_mdinfile: str,
              ll_mdinfile: str,
              target: str = ":1",
              targetcharge: int = -2):
    """
    Assume there is only one frame.
    cutoff = 4.5
    """
    mask_str = "((%s)<:%.1f) & !(%s)" % (target, cutoff, target)
    s = dpdata.System("",
        nc_file=ncfile, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
    # systems with nearby waters
    parm = load_param_file(parmfile)
    target_idx = pick_by_amber_mask(parm, target, s['coords'][0])
    nearby_idx = pick_by_amber_mask(parm, mask_str, s['coords'][0])
    parm = reorder(parm, target_idx, nearby_idx)
    # Here we need to generate 4 new parms:
    # 1. reordered, QM water
    # 2. reordered, QM water, no target
    # 3. reordered, MM water
    # 4. reordered, MM water, no target
    # parm 1
    new_nearby_idx = range(len(target_idx), len(target_idx) + len(nearby_idx))
    parm_qm = copy(parm)
    parm_qm, nwater = remove_epw(parm_qm, new_nearby_idx)
    water_mask = "(:2-%d)" % (nwater+1)
    qmmask = "(%s)|(%s)" % (target, water_mask)
    sort_parm(parm_qm)
    write_parm(parm_qm, "qmwater", hl_mdinfile, qmmask, target, charge=targetcharge)

    # parm2
    parm_qm_not = strip_atoms(parm_qm, target)
    water_mask = "(:1-%d)" % (nwater)
    sort_parm(parm_qm_not)
    write_parm(parm_qm_not, "qmwater_not", hl_mdinfile, water_mask, "")

    # parm3
    parm_mm = copy(parm)
    sort_parm(parm_mm)
    write_parm(parm_mm, "mmwater", ll_mdinfile, target, target, charge=targetcharge)

    # parm4
    parm_mm_not = strip_atoms(parm_mm, target)
    sort_parm(parm_mm_not)
    write_parm(parm_mm_not, "mmwater_not", ll_mdinfile, "", "")

    # also try to remove all epw
    parm_noep = copy(parm)
    parm_noep, _ = remove_epw(parm_noep, range(len(target_idx), s.get_natoms()))
    sort_parm(parm_mm)
    parm_noep.write_parm("noepw.parm7")


def label_atoms(parm, target_idx, qmwater_idx):
    """set qwt"""
    for ii in range(len(parm.atoms)):
        atom = parm.atoms[ii]
        if ii in target_idx:
            atom._sorted = 0
        elif ii in qmwater_idx:
            atom._sorted = 1
        else:
            atom._sorted = 2
        

def remove_epw(parm, idx):
    """Remove all EP atoms in water in idx"""
    nwater = 0
    ep = []
    for ii in idx:
        atom = parm.atoms[ii]
        res = atom.residue
        if atom.type == "EP":
            ep.append(ii)
            nwater += 1
            if res.atoms[0].name == "O":
                res.name = "qwt"
                res.atoms[0].charge += atom.charge
                #print(atom.charge, res.atoms[0].charge)
                atom.charge = 0
                # remove from res
                res.atoms.pop(res.atoms.index(atom))
            else:
                raise RuntimeError()
            # remove bonds
            for bond in atom.bonds:
                bond.delete()
            for _ in range(len(atom._bond_partners)):
                atom._bond_partners.pop()
            for _ in range(len(atom._exclusion_partners)):
                atom._exclusion_partners.pop()
        else:
            for p in atom.bonds:
                if p.atom1.type == "EP" or p.atom2.type == "EP":
                    p.delete()
            for _ in range(len(atom._exclusion_partners)):
                atom._exclusion_partners.pop()
            for p in atom.children:
                if p.type == "EP":
                    atom.children.pop(atom.children.index(p))

    mask = "@" + ",".join([ str(ii+1) for ii in ep])
    parm = strip_atoms(parm, mask)
    return parm, nwater

def strip_atoms(parm, mask):
    parm.strip(mask)
    return parm

def reorder(parm, target_idx, nearby_idx):
    label_atoms(parm, target_idx, qmwater_idx=nearby_idx)
    parm.atoms.sort(key=lambda a: a._sorted)
    return parm

def sort_parm(parm):
    # sort residues
    parm.residues.sort(key=sort_residues)

def write_parm(parm, fn, mdinfile, qmmask, target, charge=0):
    parm.write_parm("%s.parm7" % fn)
    parm.write_rst7("%s.rst7" % fn)
    PT.writeCoordinates(parm, '%s.nc' % fn).execute()
    ifqnt = str(int(target!=""))
    with open(mdinfile) as f, open("%s.mdin" % fn, 'w') as fw:
        fw.write(f.read().replace("%QMMASK%", qmmask).replace("%NOSHAKEMASK%", target).replace("%CHARGE%", str(charge)).replace("%IFQNT%", ifqnt))


def sort_residues(res):
    keys = {
        "qwt": 1,
        "WAT": 2,
    }
    return keys.get(res.name, 0)

def run(args):
    qmbuffer(cutoff=args.cutoff,
            parmfile=args.parm7_file,
            ncfile=args.nc,
            ll_mdinfile=args.ll,
            hl_mdinfile=args.hl,
            target=args.qm_region,
            targetcharge=args.charge,
            )
