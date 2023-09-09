import re

import numpy as np
import parmed.geometry


def cpt_dist_and_grd(ca, cb):
    """Return the distance between two points and the gradient.

    Parameters
    ----------
    ca : numpy.array, shape=(3,)
        A point in space

    cb : numpy.array, shape=(3,)
        A point in space

    Returns
    -------
    q : float
        The distance |a-b|

    dqda : numpy.array, shape=(3,)
        Gradient with respect to a

    dqdb : numpy.array, shape=(3,)
        Gradient with respect to b
    """
    rabv = ca - cb
    rab2 = np.dot(rabv, rabv)
    rab = np.sqrt(rab2)
    z = rab
    dzdra = (0.5 / rab) * (2 * rabv[:])
    dzdrb = (0.5 / rab) * (2 * rabv[:]) * (-1.0)
    return z, dzdra, dzdrb


def cpt_angle_and_grd(ca, cb, cc):
    """Return the angle of 3 points in degrees and gradient.

    Parameters
    ----------
    ca : numpy.array, shape=(3,)
        A point in space

    cb : numpy.array, shape=(3,)
        A point in space

    cc : numpy.array, shape=(3,)
        A point in space

    Returns
    -------
    q : float
        The angle in degrees

    dqda : numpy.array, shape=(3,)
        Gradient with respect to a

    dqdb : numpy.array, shape=(3,)
        Gradient with respect to b

    dqdc : numpy.array, shape=(3,)
        Gradient with respect to c
    """
    f = 180.0 / np.pi
    x1, y1, z1 = ca[0], ca[1], ca[2]
    x2, y2, z2 = cb[0], cb[1], cb[2]
    x3, y3, z3 = cc[0], cc[1], cc[2]
    v1 = np.array([x2 - x1, y2 - y1, z2 - z1])
    v2 = np.array([x2 - x3, y2 - y3, z2 - z3])
    l1 = np.sqrt(np.dot(v1, v1))
    l2 = np.sqrt(np.dot(v2, v2))
    n = np.dot(v1, v2)
    d = l1 * l2
    cosa = n / d
    z = np.arccos(cosa) * f
    dzdca = -f / np.sqrt(1.0 - cosa**2)
    v1 / l1
    dcadv1 = v2 / d + (-n / (d * d)) * l2 * (v1 / l1)
    dcadv2 = v1 / d + (-n / (d * d)) * l1 * (v2 / l2)

    dv1da = np.zeros((3, 3))
    dv1db = np.zeros((3, 3))

    dv1da[0, 0] = -1.0
    dv1da[1, 1] = -1.0
    dv1da[2, 2] = -1.0
    dv1db[0, 0] = 1.0
    dv1db[1, 1] = 1.0
    dv1db[2, 2] = 1.0
    dv2db = dv1db
    dv2dc = dv1da

    dzda = dzdca * np.dot(dcadv1, dv1da)
    dzdb = dzdca * (np.dot(dcadv1, dv1db) + np.dot(dcadv2, dv2db))
    dzdc = dzdca * np.dot(dcadv2, dv2dc)

    return z, dzda, dzdb, dzdc


def cpt_r12_and_grd(rstwt1, rstwt2, ca, cb, cc, cd):
    """Return a linear combination of 2 distances and gradient.

    Parameters
    ----------
    rstwt1 : float
        The weight 1

    rstwt2 : float
        The weight 2
    ca : numpy.array, shape=(3,)
        A point in space

    cb : numpy.array, shape=(3,)
        A point in space

    cc : numpy.array, shape=(3,)
        A point in space

    cd : numpy.array, shape=(3,)
        A point in space



    Returns
    -------
    q : float
        The linear combination of distances rstwt1*|a-b| + rstwt2*|c-d|

    dqda : numpy.array, shape=(3,)
        Gradient with respect to a

    dqdb : numpy.array, shape=(3,)
        Gradient with respect to b

    dqdc : numpy.array, shape=(3,)
        Gradient with respect to c

    dqdd : numpy.array, shape=(3,)
        Gradient with respect to d
    """
    rabv = ca - cb
    rcdv = cc - cd
    rab2 = np.dot(rabv, rabv)
    rcd2 = np.dot(rcdv, rcdv)
    rab = np.sqrt(rab2)
    rcd = np.sqrt(rcd2)
    z = rstwt1 * rab + rstwt2 * rcd
    dzdra = rstwt1 * ((0.5 / rab) * (2 * rabv[:]))
    dzdrb = rstwt1 * ((0.5 / rab) * (2 * rabv[:]) * (-1.0))
    dzdrc = rstwt2 * ((0.5 / rcd) * (2 * rcdv[:]))
    dzdrd = rstwt2 * ((0.5 / rcd) * (2 * rcdv[:]) * (-1.0))
    return z, dzdra, dzdrb, dzdrc, dzdrd


class Restraint:
    """Stores restraint definition.

    Parameters
    ----------
    line : str
        The fortran namelist definition of the restraint.
        The namelist should begin with &rst and end with / or &end

    Attributes
    ----------
    line : str
        The fortran namelist string used to define the restraint

    angle : bool
        Indicates if the restraint is a 3-atom angle definition

    dihed : bool
        Indicates if the restraint is a 4-atom dihedral angle definition

    r12 : bool
        Indicates if the restraint is a 4-atom R12 definition

    iat : list of int
        The indexes of the atoms involved in the restraint (1-based indexing)

    rstwt : list of float
        The coefficients of the linear combination of distances used in the
        R12 distance definition

    r1 : float
        The lower bound of the restraint, below which the penalty is linear

    r2 : float or str
        The lower center of the flat harmonic. If string, then this is a
        template definition

    r3 : float or str
        The upper center of the flat harmonic. If string, then this is a
        template definition

    r4 : float
        The upper bound of the restraint, above which the penalty is linear

    rk2 : float or str
        The quadratic penalty between r1 and r2
        Note that angle and dihedral restraints list this value in units of
        (kcal/mol) * radian**(-2); however, it is immediately converted to
        (kcal/mol) * degree**(-2) because the angles are assumed to be in
        degrees -- as per Amber convention.
        If string, then this is a template definition

    rk3 : float or str
        The quadratic penalty between r3 and r4
        Note that angle and dihedral restraints list this value in units of
        (kcal/mol) * radian**(-2); however, it is immediately converted to
        (kcal/mol) * degree**(-2) because the angles are assumed to be in
        degrees -- as per Amber convention.
        If string, then this is a template definition
    """

    def __init__(self, line: str):
        self.line = line
        self.angle = False
        self.dihed = False
        self.r12 = False
        self.iat = []
        self.rstwt = []
        self.r1 = -1000.0
        self.r2 = 0
        self.r3 = 0
        self.r4 = 1000.0
        self.rk2 = 0
        self.rk3 = 0

        line = line.replace(",", " ")
        cols = line.split()
        cols.remove("&rst")
        if "/" in cols:
            cols.remove("/")
        if "&end" in cols:
            cols.remove("&end")
        line = " ".join(cols)

        keyvals = re.findall(r"(.*?)=([^=]+)", line)
        for i in range(len(keyvals)):
            k, v = keyvals[i]
            keyvals[i] = (k.strip(), v.strip())
        for i in range(len(keyvals)):
            k, v = keyvals[i]
            if len(k) == 0:
                pk, pv = keyvals[i - 1]
                cols = pv.split()
                k = cols.pop()
                pv = " ".join(cols)
                keyvals[i - 1] = pk, pv
                keyvals[i] = k, v

        for k, v in keyvals:
            if k == "iat":
                self.iat = [int(x) for x in v.split()]
                self.angle = False
                if len(self.iat) == 3:
                    self.angle = True
                elif len(self.iat) == 4:
                    if "rstwt" in [k for k, v in keyvals]:
                        self.angle = False
                        self.r12 = True
                    else:
                        self.angle = True
                        self.dihed = True
            elif k == "rk2":
                try:
                    self.rk2 = float(v)
                except Exception:
                    self.rk2 = v
            elif k == "rk3":
                try:
                    self.rk3 = float(v)
                except Exception:
                    self.rk3 = v
            elif k == "r1":
                self.r1 = float(v)
            elif k == "r2":
                try:
                    self.r2 = float(v)
                except Exception:
                    self.r2 = v
            elif k == "r3":
                try:
                    self.r3 = float(v)
                except Exception:
                    self.r3 = v
            elif k == "r4":
                self.r4 = float(v)
            elif k == "rstwt":
                self.rstwt = [float(x) for x in v.split()]

        if self.angle:
            f = (np.pi / 180.0) ** 2
            if isinstance(self.rk2, float):
                self.rk2 *= f
            if isinstance(self.rk3, float):
                self.rk3 *= f

    def get_wilson_elements(self, crds, aidxs):
        """Return a column of the Wilson B matrix, the elements
        of which are dq/dx, where q is the restraint coordinate
        and x is a Cartesian coordinate.  The length of the output
        array is 3 * len(aidxs).

        Parameters
        ----------
        crds : numpy.array, shape=(nat,3)
            The coordinates of the system

        aidxs : list of int
            Each element is the 0-based index of the atom in
            the full system.  In contrast, crds may be a
            petite list of coordinates, where the petite list
            is len(aidxs)

        Returns
        -------
        q : float
            The restraint coordinate value

        B : numpy.array, shape=3*len(aidxs)
            The values of dq/dx for this restraint
        """
        q = None
        dqdx = np.zeros(3 * len(aidxs))

        umap = {}
        for u, a in enumerate(aidxs):
            umap[a] = u

        if len(self.iat) == 2:
            ia = self.iat[0] - 1
            if ia in umap:
                ia = umap[ia]
            else:
                raise Exception("Missing %i in aidxs" % (ia))

            ib = self.iat[1] - 1
            if ib in umap:
                ib = umap[ib]
            else:
                raise Exception("Missing %i in aidxs" % (ib))

            q, dqda, dqdb = cpt_dist_and_grd(crds[ia, :], crds[ib, :])

            for k in range(3):
                dqdx[k + ia * 3] = dqda[k]
                dqdx[k + ib * 3] += dqdb[k]

        elif self.dihed:
            ia = self.iat[0] - 1
            if ia in umap:
                ia = umap[ia]
            else:
                raise Exception("Missing %i in aidxs" % (ia))

            ib = self.iat[1] - 1
            if ib in umap:
                ib = umap[ib]
            else:
                raise Exception("Missing %i in aidxs" % (ib))

            ic = self.iat[2] - 1
            if ic in umap:
                ic = umap[ic]
            else:
                raise Exception("Missing %i in aidxs" % (ic))

            ie = self.iat[3] - 1
            if ie in umap:
                ie = umap[ie]
            else:
                raise Exception("Missing %i in aidxs" % (ie))

            q = parmed.geometry.dihedral(
                crds[ia, :], crds[ib, :], crds[ic, :], crds[ie, :]
            )  # * np.pi / 180.
            raise Exception("Dihedrals not implemented")

        elif self.angle:
            ia = self.iat[0] - 1
            if ia in umap:
                ia = umap[ia]
            else:
                raise Exception("Missing %i in aidxs" % (ia))

            ib = self.iat[1] - 1
            if ib in umap:
                ib = umap[ib]
            else:
                raise Exception("Missing %i in aidxs" % (ib))

            ic = self.iat[2] - 1
            if ic in umap:
                ic = umap[ic]
            else:
                raise Exception("Missing %i in aidxs" % (ic))

            q, dqda, dqdb, dqdc = cpt_angle_and_grd(
                crds[ia, :], crds[ib, :], crds[ic, :]
            )

            for k in range(3):
                dqdx[k + ia * 3] = dqda[k]
                dqdx[k + ib * 3] += dqdb[k]
                dqdx[k + ic * 3] += dqdc[k]

        elif len(self.iat) == 4:
            ia = self.iat[0] - 1
            if ia in umap:
                ia = umap[ia]
            else:
                raise Exception("Missing %i in aidxs" % (ia))

            ib = self.iat[1] - 1
            if ib in umap:
                ib = umap[ib]
            else:
                raise Exception("Missing %i in aidxs" % (ib))

            ic = self.iat[2] - 1
            if ic in umap:
                ic = umap[ic]
            else:
                raise Exception("Missing %i in aidxs" % (ic))

            ie = self.iat[3] - 1
            if ie in umap:
                ie = umap[ie]
            else:
                raise Exception("Missing %i in aidxs" % (ie))

            q, dqda, dqdb, dqdc, dqde = cpt_r12_and_grd(
                self.rstwt[0],
                self.rstwt[1],
                crds[ia, :],
                crds[ib, :],
                crds[ic, :],
                crds[ie, :],
            )

            for k in range(3):
                dqdx[k + ia * 3] = dqda[k]
                dqdx[k + ib * 3] += dqdb[k]
                dqdx[k + ic * 3] += dqdc[k]
                dqdx[k + ie * 3] += dqde[k]
        return q, dqdx


class Disang:
    """A collection of restraints.

    Parameters
    ----------
    fname : str
        The name of the Amber restraint file (disang file)

    Attributes
    ----------
    restraints : list of Restraint
        The restraint definitions
    """

    def __init__(self, fname: str) -> None:
        restraints = []
        rline = ""
        with open(fname) as fh:
            for line in fh:
                line = line.split("!")[0].strip()
                if "&rst" in line:
                    rline = line
                    if not ("/" in line or "&end" in line):
                        for line in fh:
                            line = line.split("!")[0].strip()
                            if "/" in line or "&end" in line:
                                rline += " " + line
                                break
                            else:
                                rline += " " + line
                    rline = rline.replace("/", " /")
                    restraints.append(Restraint(rline))
                    rline = ""
        self.restraints = restraints

    def get_unique_atom_idxs(self):
        """Return a list of unique atom indexes used to define
        retraints.

        Returns
        -------
        aidxs : list of int
            The sorted list of 0-based atom indexes
        """
        aidxs = []
        for res in self.restraints:
            aidxs.extend([i - 1 for i in res.iat])

        return np.unique(aidxs)

    def get_drdq(self, crds_all: np.ndarray) -> np.ndarray:
        nall = crds_all.shape[0]

        # aidxs is a petite array of N' atom indexes
        aidxs = self.get_unique_atom_idxs()
        # rcs is an array of Nrc generalized coordinates
        # B is a Nrc x 3N' matrix of derivatives
        nres = len(self.restraints)

        nat = len(aidxs)
        B = np.zeros((nres, 3 * nat))

        crds = crds_all[aidxs, :]

        qs = np.zeros((nres,))
        for ires in range(nres):
            q, col = self.restraints[ires].get_wilson_elements(crds, aidxs)
            qs[ires] = q
            B[ires, :] = col[:]

        # Binv is the 3N' x Nrc pseudoinverse matrix
        # used to project the forces
        drdq = np.linalg.pinv(B)

        drdq_all = np.zeros((nall, 3, nres))
        drdq_all[aidxs, :, :] = drdq.reshape((nat, 3, nres))
        return drdq_all
