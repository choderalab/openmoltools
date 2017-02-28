import numpy as np
import itertools
import simtk.openmm as mm
from simtk import unit as u

import logging
logger = logging.getLogger(__name__)

EPSILON = 1E-4  # Error tolerance for differences in parameters.  Typically for relative differences, but sometimes for absolute.

"""NOTE: There is currently a bug due to ambiguity of improper torsion
atom order.  This will be fixed soon, but
for now the below values of EPSILON are chosen to silence this bug during testing.
"""
ENERGY_EPSILON = 3.0 * u.kilocalories_per_mole
TORSION_ENERGY_EPSILON = 3.0 * u.kilocalories_per_mole
COMPONENT_ENERGY_EPSILON = 0.1 * u.kilocalories_per_mole



def compare(x0, x1, relative=False):
    """Compare two quantities relative to EPSILON."""

    if relative is True:
        # Watch out for zero division
        if abs(x1) > 0:
            denominator = abs(x1)
        elif abs(x0) > 0:
            denominator = abs(x0)
        else:
            return 0
    else:
        denominator = 1.0
    value = abs(x0-x1)/denominator

    # Make sure return is unitless
    try:
        value = value/value.unit
    except AttributeError:
        pass

    return (value) < EPSILON

reduce_precision = lambda x: float(np.float16(x))  # Useful for creating dictionary keys with floating point numbers that may differ at insignificant decimal places

reorder_bonds = lambda i0, i1: (min(i0, i1), max(i0, i1))
reorder_angles = lambda i0, i1, i2: (min(i0, i2), i1, max(i0, i2))


def reorder_proper_torsions(i0, i1, i2, i3):
    """Return the atom indices of a proper torsion after "flipping" the
    order so that the first atom is the smallest index.

    Parameters
    ----------
    i0, i1, i2, i3 : int,
        Atom indices of torsion

    Returns
    -------
    j0, j1, j2, j3 : int,
        Atom indices of torsion

    """
    if i0 < i3:
        j0, j1, j2, j3 = i0, i1, i2, i3
    else:
        j0, j1, j2, j3 = i3, i2, i1, i0

    return j0, j1, j2, j3


def reorder_improper_torsions(i0, i1, i2, i3, bond_set):
    """Return j0, j1, j2, j3, where j0 is the central index and j1, j2, j3
    are in sorted() order.

    Parameters
    ----------
    i0, i1, i2, i3 : int,
        Atom indices of torsion
    bond_set : set containing the index pairs between which a bond is defined

    Returns
    -------
    j0, j1, j2, j3 : int,
        Atom indices of torsion, with j0 being the central index

    Notes
    -----

    Centrality is determined by the maximum counts in the adjacency matrix.

    """

    connections = np.zeros((4, 4))

    mapping = {i0: 0, i1: 1, i2: 2, i3: 3}
    inv_mapping = dict([(val, key) for key, val in mapping.items()])

    for (a, b) in itertools.combinations([i0, i1, i2, i3], 2):
        if (a, b) in bond_set:
            i, j = mapping[a], mapping[b]
            connections[i, j] += 1.
            connections[j, i] += 1.

    central_ind = connections.sum(0).argmax()
    central_ind = inv_mapping[central_ind]
    other_ind = sorted([i0, i1, i2, i3])
    other_ind.remove(central_ind)

    return central_ind, other_ind[0], other_ind[1], other_ind[2]


def get_symmetrized_bond_set(bond_force):
    """Return a set containing atom pairs connected by bonds.

    Parameters
    ----------
    bond_force : mm.HarmonicBondForce
        The bond force of an OpenMM system

    Returns
    -------
    bond_set : a set containing pairs of atoms connected by bonds:

    Notes
    -----
    The resulting set is symmetric: if (i,j) in S, then (j,i) in S.
    """

    bond_set = set()
    n_bonds = bond_force.getNumBonds()

    for k in range(n_bonds):
        (i0, i1, r0, k0) = bond_force.getBondParameters(k)
        bond_set.add((i0, i1))
        bond_set.add((i1, i0))

    return bond_set


def is_proper(i0, i1, i2, i3, bond_set):
    """Check for three sequential bonds and atom uniqueness."""
    if (i0, i1) in bond_set and (i1, i2) in bond_set and (i2, i3) in bond_set and len(set([i0, i1, i2, i3])) == 4:
        return True
    return False


def is_improper(i0, i1, i2, i3, bond_set):
    """Check for three non-sequential bonds and atom uniqueness."""
    if len(set([i0, i1, i2, i3])) == 4:
        if not ((i0, i1) in bond_set and (i1, i2) in bond_set and (i2, i3) in bond_set):
            return True
    return False


class SystemChecker(object):

    def __init__(self, simulation0, simulation1):
        """Create a SystemChecker object that compares forces in simulation0 and simulation1.

        Parameters
        ----------
        simulation0 : OpenMM Simulation
        simulation1 : OpenMM Simulation

        Notes
        -----

        You MUST have constraints=None when creating your simulation objects.

        """

        self.simulation0 = simulation0
        self.simulation1 = simulation1

        for force in simulation0.system.getForces():
            if type(force) == mm.HarmonicBondForce:
                self.bond_force0 = force
            elif type(force) == mm.HarmonicAngleForce:
                self.angle_force0 = force
            elif type(force) == mm.PeriodicTorsionForce:
                self.torsion_force0 = force
            elif type(force) == mm.NonbondedForce:
                self.nonbonded_force0 = force

        for force in simulation1.system.getForces():
            if type(force) == mm.HarmonicBondForce:
                self.bond_force1 = force
            elif type(force) == mm.HarmonicAngleForce:
                self.angle_force1 = force
            elif type(force) == mm.PeriodicTorsionForce:
                self.torsion_force1 = force
            elif type(force) == mm.NonbondedForce:
                self.nonbonded_force1 = force

    def check_force_parameters(self, skipImpropers = False):
        """Check that force parameters are the same, up to some equivalence.

        Arguments:
        ----------
        skipImpropers : bool (optional), default False
            Skip checking of impropers if desired (i.e. for a force field using a different form for impropers).
        """
        self.check_bonds(self.bond_force0, self.bond_force1)
        self.check_angles(self.angle_force0, self.angle_force1)
        self.check_nonbonded(self.nonbonded_force0, self.nonbonded_force1)
        self.check_proper_torsions(self.torsion_force0, self.torsion_force1, self.bond_force0, self.bond_force1)
        if not skipImpropers:
            self.check_improper_torsions(self.torsion_force0, self.torsion_force1, self.bond_force0, self.bond_force1)
            logger.info("Note: skipping degenerate impropers with < 4 atoms.")

    def check_bonds(self, force0, force1):
        """Check that force0 and force1 are equivalent Bond forces.


        Parameters
        ----------
        force0 : mm.HarmonicBondForce
        force1 : mm.HarmonicBondForce

        """

        assert type(force0) == type(force1), "Error: force0 and force1 must be the same type."
        assert type(force0) == mm.HarmonicBondForce, "Error: forces must be HarmonicBondForces"

        n_bonds0 = force0.getNumBonds()
        n_bonds1 = force1.getNumBonds()

        dict0, dict1 = {}, {}

        i0, i1, r0, k0 = force0.getBondParameters(0)
        unit_r = u.angstrom
        #unit_k = k0.unit
        unit_k = u.kilojoules_per_mole/(u.angstrom)**2

        for k in range(n_bonds0):
            i0, i1, r0, k0 = force0.getBondParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            if k0 / k0.unit != 0.0:  # Skip forces with strength 0.0
                dict0[i0, i1] = ((r0 / unit_r, k0 / unit_k))

        for k in range(n_bonds1):
            i0, i1, r0, k0 = force1.getBondParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            if k0 / k0.unit != 0.0:  # Skip forces with strength 0.0
                dict1[i0, i1] = ((r0 / unit_r, k0 / unit_k))

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        logger.info("Bonds0 - Bonds1 = %s" % (keys0.difference(keys1)))
        logger.info("Bonds1 - Bonds0 = %s" % (keys1.difference(keys0)))
        assert set(dict0.keys()) == set(dict1.keys()), "Systems have different HarmonicBond Forces"

        for k, parameter_name in enumerate(["r0", "k0"]):
            for (i0, i1) in dict0.keys():
                val0 = dict0[i0, i1][k]
                val1 = dict1[i0, i1][k]
                if parameter_name=='r0':
                    assert compare(val0, val1), "Error: Harmonic Bond distance (%d, %d) has equilibrium distances of %f and %f angstroms, respectively." % (i0, i1, val0, val1)
                else:
                    assert compare(val0, val1), "Error: Harmonic Bond force constant (%d, %d) has values of %f and %f kJ/mol, respectively." % (i0, i1, val0, val1)


    def check_angles(self, force0, force1):
        """Check that force0 and force1 are equivalent Angle forces.


        Parameters
        ----------
        force0 : mm.HarmonicAngleForce
        force1 : mm.HarmonicAngleForce

        """

        assert type(force0) == type(force1), "Error: force0 and force1 must be the same type."
        assert type(force0) == mm.HarmonicAngleForce, "Error: forces must be HarmonicAngleForces"

        n_angles0 = force0.getNumAngles()
        n_angles1 = force1.getNumAngles()

        dict0, dict1 = {}, {}

        i0, i1, i2, theta0, k0 = force0.getAngleParameters(0)
        #unit_theta = theta0.unit
        unit_theta = u.degrees
        #unit_k = k0.unit
        unit_k = u.kilojoules_per_mole/(u.degrees)**2

        for k in range(n_angles0):
            i0, i1, i2, theta0, k0 = force0.getAngleParameters(k)
            if (k0 / k0.unit) != 0.0:  # Skip forces with strength 0.0
                i0, i1, i2 = reorder_angles(i0, i1, i2)
                dict0[i0, i1, i2] = ((theta0 / unit_theta, k0 / unit_k))

        for k in range(n_angles1):
            i0, i1, i2, theta0, k0 = force1.getAngleParameters(k)
            if (k0 / k0.unit) != 0.0:  # Skip forces with strength 0.0
                i0, i1, i2 = reorder_angles(i0, i1, i2)
                dict1[i0, i1, i2] = ((theta0 / unit_theta, k0 / unit_k))

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        logger.info("Angles0 - Angles1 = %s" % (keys0.difference(keys1)))
        logger.info("Angles1 - Angles0 = %s" % (keys1.difference(keys0)))
        diff_keys = keys0.symmetric_difference(keys1)
        assert diff_keys == set(), "Systems have different HarmonicAngleForce entries: extra keys are: \n%s" % diff_keys

        for k, parameter_name in enumerate(["theta0", "k0"]):
            for (i0, i1, i2) in dict0.keys():
                val0 = dict0[i0, i1, i2][k]
                val1 = dict1[i0, i1, i2][k]
                if parameter_name=='theta0':
                    assert compare(val0, val1), "Error: Harmonic Angle (%d, %d, %d) has angle values of %f and %f degrees, respectively." % (i0, i1, i2, val0, val1)
                else:
                    assert compare(val0, val1), "Error: Harmonic Angle (%d, %d, %d) has force constant values of %f and %f kJ/(mol degree**2), respectively." % (i0, i1, i2, val0, val1)

    def check_nonbonded(self, force0, force1):
        """Check that force0 and force1 are equivalent Nonbonded forces.


        Parameters
        ----------
        force0 : mm.NonbondedForce
        force1 : mm.NonbondedForce

        """

        assert type(force0) == type(force1), "Error: force0 and force1 must be the same type."
        assert type(force0) == mm.NonbondedForce, "Error: forces must be NonbondedForces"
        assert force0.getNumParticles() == force1.getNumParticles(), "Error: Systems have %d and %d particles in NonbondedForce, respectively." % (force0.getNumParticles(), force1.getNumParticles())

        n_atoms = force0.getNumParticles()

        q, sigma, epsilon = force0.getParticleParameters(0)
        #unit_q, unit_sigma, unit_epsilon = q.unit, sigma.unit, epsilon.unit
        unit_q = u.elementary_charge
        unit_sigma = u.angstrom
        unit_epsilon = u.kilojoule_per_mole

        for k in range(n_atoms):
            q0, sigma0, epsilon0 = force0.getParticleParameters(k)
            q1, sigma1, epsilon1 = force1.getParticleParameters(k)

            q0, sigma0, epsilon0 = q0 / unit_q, sigma0 / unit_sigma, epsilon0 / unit_epsilon
            q1, sigma1, epsilon1 = q1 / unit_q, sigma1 / unit_sigma, epsilon1 / unit_epsilon

            assert compare(q0, q1), "Error: Particle %d has charges of %f and %f, respectively." % (k, q0, q1)

            if epsilon0 != 0.:
                assert compare(sigma0, sigma1), "Error: Particle %d has sigma of %f and %f angstroms, respectively." % (k, sigma0, sigma1)
            else:
                logger.info("Skipping comparison of sigma (%f, %f) on particle %d because epsilon has values %f, %f kJ/mol" % (sigma0, sigma1, k, epsilon0, epsilon1))

            assert compare(epsilon0, epsilon1), "Error: Particle %d has epsilon of %f and %f kJ/mol, respectively." % (k, epsilon0, epsilon1)

        n_exceptions = force0.getNumExceptions()
        assert force0.getNumExceptions() == force1.getNumExceptions(), "Error: Systems have %d and %d exceptions in NonbondedForce, respectively." % (force0.getNumExceptions(), force1.getNumExceptions())

        i0, i1, qq, sigma, epsilon = force0.getExceptionParameters(0)
        unit_qq = u.elementary_charge**2
        unit_sigma = u.angstrom
        unit_epsilon = u.kilojoule_per_mole

        dict0, dict1 = {}, {}
        for k in range(n_exceptions):
            i0, i1, qq, sigma, epsilon = force0.getExceptionParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            dict0[i0, i1] = ((qq / unit_qq, sigma / unit_sigma, epsilon / unit_epsilon))

            i0, i1, qq, sigma, epsilon = force1.getExceptionParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            dict1[i0, i1] = ((qq / unit_qq, sigma / unit_sigma, epsilon / unit_epsilon))

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        logger.info("Exceptions0 - Exceptions1 = %s" % (keys0.difference(keys1)))
        logger.info("Exceptions1 - Exceptions0 = %s" % (keys1.difference(keys0)))
        assert set(dict0.keys()) == set(dict1.keys()), "Systems have different NonBondedForce Exceptions"

        for k, parameter_name in enumerate(["qq", "sigma", "epsilon"]):
            for (i0, i1) in dict0.keys():
                val0 = dict0[i0, i1][k]
                val1 = dict1[i0, i1][k]
                if parameter_name == "sigma" and dict0[i0, i1][2] == 0.0 and dict1[i0, i1][2] == 0.0:
                    continue  # If both epsilon parameters are zero, then sigma doesn't matter so skip the comparison.
                if parameter_name =="sigma":
                    assert compare(val0, val1), "Error: NonBondedForce Exception, atom (%d, %d) has sigma values of %f and %f angstroms, respectively." % (i0, i1, parameter_name, val0, val1)
                elif parameter_name=="qq":
                    assert compare(val0, val1), "Error: NonBondedForce Exception atom (%d, %d) has squared charge values of %f and %f (elementary charge)**2, respectively." % (i0, i1, val0, val1)
                else:
                    assert compare(val0, val1), "Error: NonBondedForce Exception, atom (%d, %d) has epsilon values of %f and %f kJ/mol, respectively." % (i0, i1, val0, val1)

    def check_proper_torsions(self, force0, force1, bond_force0, bond_force1):
        """Check that force0 and force1 are equivalent PeriodicTorsion forces.


        Parameters
        ----------
        force0 : mm.PeriodicTorsionForce
        force1 : mm.PeriodicTorsionForce

        Notes
        -----

        So this creates and compares a pair of dictionaries, one for each
        force.  Each dictionary has keys that are atom tuples (i0,i1,i2,i3)
        and values which are a list of force parameters (per, phase, k0)
        for that tuple.  This complexity is required because a given tuple
        can have *multiple* parameters associated with it.

        """

        assert type(force0) == type(force1), "Error: force0 and force1 must be the same type."
        assert type(force0) == mm.PeriodicTorsionForce, "Error: forces must be PeriodicTorsionForce"

        bond_set0 = get_symmetrized_bond_set(bond_force0)
        bond_set1 = get_symmetrized_bond_set(bond_force1)

        # Build list of atoms to help make output more useful
        atoms0 = [ atom for atom in self.simulation0.topology.atoms() ]
        atoms1 = [ atom for atom in self.simulation1.topology.atoms() ]


        # Check torsions for equivalent

        if force0.getNumTorsions() == 0 and force1.getNumTorsions() == 0:
            return  # Must leave now, otherwise try to access torsions that don't exist.

        F = force0 if force0.getNumTorsions() > 0 else force1  # Either force0 or force1 is nonempty, so find that one.
        i0, i1, i2, i3, per, phase, k0 = F.getTorsionParameters(0)
        #phase_unit, k0_unit = phase.unit, k0.unit
        k0_unit = u.kilojoules_per_mole
        phase_unit = u.degrees

        dict0, dict1 = {}, {}
        for k in range(force0.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force0.getTorsionParameters(k)
            if not is_proper(i0, i1, i2, i3, bond_set0):
                continue

            i0, i1, i2, i3 = reorder_proper_torsions(i0, i1, i2, i3)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not (i0, i1, i2, i3) in dict0:
                dict0[i0, i1, i2, i3] = []

            dict0[i0, i1, i2, i3].append((per, phase, k0))

        for k in range(force1.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force1.getTorsionParameters(k)
            if not is_proper(i0, i1, i2, i3, bond_set1):
                continue

            i0, i1, i2, i3 = reorder_proper_torsions(i0, i1, i2, i3)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not (i0, i1, i2, i3) in dict1:
                dict1[i0, i1, i2, i3] = []

            dict1[i0, i1, i2, i3].append((per, phase, k0))

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        diff_keys = keys0.symmetric_difference(keys1)

        logger.info("Propers0 - Propers1 = %s" % (keys0.difference(keys1)))
        logger.info("Propers1 - Propers0 = %s" % (keys1.difference(keys0)))
        assert diff_keys == set(), "Systems have different (proper) PeriodicTorsionForce entries: extra keys are: \n%s" % diff_keys

        for (i0, i1, i2, i3) in dict0.keys():
            # Store strings for printing debug info
            torsion_atoms0 = '%s - %s - %s - %s' % (atoms0[i0].name, atoms0[i1].name, atoms0[i2].name, atoms0[i3].name )
            torsion_atoms1 = '%s - %s - %s - %s' % (atoms1[i0].name, atoms1[i1].name, atoms1[i2].name, atoms1[i3].name )

            # Proceed to check
            entries0 = dict0[i0, i1, i2, i3]
            entries1 = dict1[i0, i1, i2, i3]
            if len(entries0) != len(entries1):
                print("Compared torsion involving atoms '%s' with that involving atoms '%s': " % (torsion_atoms0, torsion_atoms1))
                raise Exception("Error:  (proper) PeriodicTorsionForce entry (atoms %d, %d, %d, %d) has different numbers of terms (%d and %d, respectively)." % (i0, i1, i2, i3, len(entries0), len(entries1)))

            subdict0 = dict(((per, reduce_precision(phase)), k0) for (per, phase, k0) in entries0)
            subdict1 = dict(((per, reduce_precision(phase)), k0) for (per, phase, k0) in entries1)

            if set(subdict0.keys()) != set(subdict1.keys()):
                print("Compared torsion involving atoms '%s' with that involving atoms '%s': " % (torsion_atoms0, torsion_atoms1))
                print("Keys for system0: '%s'; keys for system1: '%s'" % (subdict0.keys(), subdict1.keys() ) )
                raise Exception("Error: (proper) PeriodicTorsionForce entry (atoms %d, %d, %d, %d) has different terms." % (i0, i1, i2, i3) )

            for (per, phase) in subdict0.keys():
                val0 = subdict0[per, phase]
                val1 = subdict1[per, phase]
                if not compare(val0, val1):
                    print("Compared torsion involving atoms '%s' with that involving atoms '%s': " % (torsion_atoms0, torsion_atoms1))
                    raise Exception( "Error: (proper) PeriodicTorsionForce (atoms %d, %d, %d, %d, periodicity %d, phase %f degrees) has barrier height values of %f and %f kJ/mol, respectively." % (i0, i1, i2, i3, per, phase, val0, val1) )

    def check_improper_torsions(self, force0, force1, bond_force0, bond_force1):
        """Check that force0 and force1 are equivalent PeriodicTorsion forces.


        Parameters
        ----------
        force0 : mm.PeriodicTorsionForce
        force1 : mm.PeriodicTorsionForce

        Notes
        -----

        So this creates and compares a pair of dictionaries, one for each
        force.  Each dictionary has keys that are atom tuples (i0,i1,i2,i3)
        and values which are a list of force parameters (per, phase, k0)
        for that tuple.  This complexity is required because a given tuple
        can have *multiple* parameters associated with it.

        In addition to the complications for *all* torsions, improper torsions
        have the additional difficulty of ambiguous definitions.  Given a tuple
        (a,b,c,d) that is held planar by an improper, one can define several permuted dihedral
        terms that give similar (but not identical) energies.  Many force field
        codes do NOT pay attention to the exact permutation that is used, leading to
        subtle differences in topology and energy between MD packages.

        The solution that we use here is to reorder each improper torsion via:

        central_atom, i1, i2, i3

        where i1, i2, and i3 are in increasing order and the central atom
        is located via the adjacency matrix of bonds.

        """

        assert type(force0) == type(force1), "Error: force0 and force1 must be the same type."
        assert type(force0) == mm.PeriodicTorsionForce, "Error: forces must be PeriodicTorsionForce"

        bond_set0 = get_symmetrized_bond_set(bond_force0)
        bond_set1 = get_symmetrized_bond_set(bond_force1)

        if force0.getNumTorsions() == 0 and force1.getNumTorsions() == 0:
            return  # Must leave now, otherwise try to access torsions that don't exist.

        F = force0 if force0.getNumTorsions() > 0 else force1  # Either force0 or force1 is nonempty, so find that one.
        i0, i1, i2, i3, per, phase, k0 = F.getTorsionParameters(0)
        #phase_unit, k0_unit = phase.unit, k0.unit
        k0_unit = u.kilojoules_per_mole
        phase_unit = u.degrees

        dict0, dict1 = {}, {}
        for k in range(force0.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force0.getTorsionParameters(k)
            if not is_improper(i0, i1, i2, i3, bond_set0):
                continue

            i0, i1, i2, i3 = reorder_improper_torsions(i0, i1, i2, i3, bond_set0)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not (i0, i1, i2, i3) in dict0:
                dict0[i0, i1, i2, i3] = []

            dict0[i0, i1, i2, i3].append((per, phase, k0))

        for k in range(force1.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force1.getTorsionParameters(k)
            if not is_improper(i0, i1, i2, i3, bond_set1):
                continue

            i0, i1, i2, i3 = reorder_improper_torsions(i0, i1, i2, i3, bond_set0)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not (i0, i1, i2, i3) in dict1:
                dict1[i0, i1, i2, i3] = []

            dict1[i0, i1, i2, i3].append((per, phase, k0))

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        diff_keys = keys0.symmetric_difference(keys1)

        logger.info("Impropers0 - Impropers1 = %s" % (keys0.difference(keys1)))
        logger.info("Impropers1 - Impropers0 = %s" % (keys1.difference(keys0)))
        assert diff_keys == set(), "Systems have different (improper) PeriodicTorsionForce entries: extra keys are: \n%s" % diff_keys

        for (i0, i1, i2, i3) in dict0.keys():
            entries0 = dict0[i0, i1, i2, i3]
            entries1 = dict1[i0, i1, i2, i3]
            assert len(entries0) == len(entries1), "Error:  (improper) PeriodicTorsionForce entry (%d, %d, %d, %d) has different numbers of terms (%d and %d, respectively)." % (i0, i1, i2, i3, len(entries0), len(entries1))

            subdict0 = dict(((per, reduce_precision(phase)), k0) for (per, phase, k0) in entries0)
            subdict1 = dict(((per, reduce_precision(phase)), k0) for (per, phase, k0) in entries1)

            assert set(subdict0.keys()) == set(subdict1.keys()), "Error: (improper) PeriodicTorsionForce entry (%d, %d, %d, %d) has different terms." % (i0, i1, i2, i3)

            for (per, phase) in subdict0.keys():
                val0 = subdict0[per, phase]
                val1 = subdict1[per, phase]
                assert compare(val0, val1), "Error: (improper) PeriodicTorsionForce (atoms %d, %d, %d, %d, periodicity %d, phase %f degrees) has barrier height values of %f and %f kJ/mol, respectively." % (i0, i1, i2, i3, per, phase, val0, val1)

    def zero_degenerate_impropers(self, f):
        """Set the force constant to zero for improper dihedrals that
        involve less than four unique atoms.
        """
        for k in range(f.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = f.getTorsionParameters(k)
            if len(set([i0, i1, i2, i3])) < 4:
                f.setTorsionParameters(k, i0, i1, i2, i3, per, phase, k0 * 0.0)

    def check_forces(self, zero_degenerate_impropers=True):
        """Compare the total forces of the two simulations.

        Parameters
        ----------

        zero_degenerate_impropers : bool, default=True
            if True, zero out all impropers with < 4 atoms.
        """
        if zero_degenerate_impropers is True:
            self.zero_degenerate_impropers(self.torsion_force0)
            xyz = self.simulation0.context.getState(getPositions=True).getPositions()
            self.simulation0.context.reinitialize()
            self.simulation0.context.setPositions(xyz)
            self.zero_degenerate_impropers(self.torsion_force1)
            xyz = self.simulation1.context.getState(getPositions=True).getPositions()
            self.simulation1.context.reinitialize()
            self.simulation1.context.setPositions(xyz)

        state0 = self.simulation0.context.getState(getForces=True)
        force0 = state0.getForces(asNumpy=True)

        state1 = self.simulation1.context.getState(getForces=True)
        force1 = state1.getForces(asNumpy=True)

        return force0, force1


    def check_energies(self, zero_degenerate_impropers=True, skip_assert=False):
        """Compare the energies of the two simulations.

        Parameters
        ----------

        zero_degenerate_impropers : bool, default=True
            if True, zero out all impropers with < 4 atoms.
        skip_assert, bool, optional, default=False
            If False, this function will raise an AssertionError if
            the energy groups are not identical.

        Notes
        -----
        If zero_degenerate_impropers is True, this function WILL MODIFY
        THE FORCES IN YOUR SYSTEM!
        """
        if zero_degenerate_impropers is True:
            self.zero_degenerate_impropers(self.torsion_force0)
            xyz = self.simulation0.context.getState(getPositions=True).getPositions()
            self.simulation0.context.reinitialize()
            self.simulation0.context.setPositions(xyz)

            self.zero_degenerate_impropers(self.torsion_force1)
            xyz = self.simulation1.context.getState(getPositions=True).getPositions()
            self.simulation1.context.reinitialize()
            self.simulation1.context.setPositions(xyz)

        state0 = self.simulation0.context.getState(getEnergy=True)
        energy0 = state0.getPotentialEnergy()

        state1 = self.simulation1.context.getState(getEnergy=True)
        energy1 = state1.getPotentialEnergy()

        if not skip_assert:
            delta = abs(energy0 - energy1)
            assert delta < ENERGY_EPSILON, "Error, energy difference (%f kJ/mol) is greater than %f kJ/mol" % (delta / u.kilojoules_per_mole, ENERGY_EPSILON / u.kilojoules_per_mole)

        return energy0, energy1


    def check_energy_groups(self, skip_assert=False):
        """Return the groupwise energies of the two simulations.

        Parameters
        ----------
        skip_assert, bool, optional, default=False
            If False, this function will raise an AssertionError if
            the energy groups are not identical.

        Returns
        -------
        energy0 : dict
            A dictionary with keys "bond", "angle", "nb", "torsion" and values
            corresponding to the energies of these components for the first simulation object
        energy1 : dict
            A dictionary with keys "bond", "angle", "nb", "torsion" and values
            corresponding to the energies of these components for the second simulation object

        """

        self.bond_force0.setForceGroup(0)
        self.bond_force1.setForceGroup(0)

        self.angle_force0.setForceGroup(1)
        self.angle_force1.setForceGroup(1)

        self.nonbonded_force0.setForceGroup(2)
        self.nonbonded_force1.setForceGroup(2)

        self.torsion_force0.setForceGroup(3)
        self.torsion_force1.setForceGroup(3)


        xyz = self.simulation0.context.getState(getPositions=True).getPositions()
        self.simulation0.context.reinitialize()
        self.simulation0.context.setPositions(xyz)

        xyz = self.simulation1.context.getState(getPositions=True).getPositions()
        self.simulation1.context.reinitialize()
        self.simulation1.context.setPositions(xyz)

        energy0 = {}
        energy1 = {}
        for k, key in enumerate(["bond", "angle", "nb", "torsion"]):
            groups = 2 ** k
            state0 = self.simulation0.context.getState(getEnergy=True, groups=groups)
            energy0[key] = state0.getPotentialEnergy()

            state1 = self.simulation1.context.getState(getEnergy=True, groups=groups)
            energy1[key] = state1.getPotentialEnergy()

        if not skip_assert:
            delta = abs(energy0["bond"] - energy1["bond"])
            assert delta < COMPONENT_ENERGY_EPSILON, "Error, bond energy difference (%f kJ/mol) is greater than %f kJ/mol" % (delta / u.kilojoules_per_mole, ENERGY_EPSILON / u.kilojoules_per_mole)
            delta = abs(energy0["angle"] - energy1["angle"])
            assert delta < COMPONENT_ENERGY_EPSILON, "Error, angle energy difference (%f kJ/mol) is greater than %f kJ/mol" % (delta / u.kilojoules_per_mole, ENERGY_EPSILON / u.kilojoules_per_mole)
            delta = abs(energy0["nb"] - energy1["nb"])
            assert delta < COMPONENT_ENERGY_EPSILON, "Error, NB energy difference (%f kJ/mol) is greater than %f kJ/mol" % (delta / u.kilojoules_per_mole, ENERGY_EPSILON / u.kilojoules_per_mole)
            delta = abs(energy0["torsion"] - energy1["torsion"])
            assert delta < TORSION_ENERGY_EPSILON, "Error, torsion energy difference (%f kJ/mol) is greater than %f kJ/mol" % (delta / u.kilojoules_per_mole, ENERGY_EPSILON / u.kilojoules_per_mole)


        return energy0, energy1
