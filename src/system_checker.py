import simtk.unit as u
import numpy as np

reduce_precision = lambda x: float(np.float16(x))  # Useful for creating dictionary keys with floating point numbers that may differ at insignificant decimal places

reorder_bonds = lambda i0, i1: (min(i0, i1), max(i0, i1))
reorder_angles = lambda i0, i1, i2: (min(i0, i2), i1, max(i0, i2))
#reorder_torsions = lambda i0, i1, i2, i3: (i0, i1, i2, i3)

def reorder_proper_torsions(i0, i1, i2, i3):
    if i0 < i3:
        j0, j1, j2, j3 = i0, i1, i2, i3
    else:
        j0, j1, j2, j3 = i3, i2, i1, i0        

    return j0, j1, j2, j3

def reorder_improper_torsions(i0, i1, i2, i3):
    """Assumes that i2 is the central atom!!!"""
    j2 = i2
    j0, j1, j3 = sorted((i0, i1, i3))
    return j0, j1, j2, j3

EPSILON = 1E-4

def get_symmetrized_bond_set(bond_force):
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

class SystemChecker(object):
    def __init__(self, system0, system1):
        self.forces0 = system0.getForces()
        self.forces1 = system1.getForces()
        self.forces0.sort(key=lambda x: type(x))
        self.forces1.sort(key=lambda x: type(x))
        
    def check_forces(self):
        self.check_bonds()
        self.check_angles()
        self.check_nonbonded()
        self.check_proper_torsions()
        self.check_improper_torsions()

    def check_bonds(self):
        force0 = self.forces0[4]  # Fix hardcoded lookup later.
        force1 = self.forces1[4]
    
        assert force0.getNumBonds() == force1.getNumBonds(), "Error: Systems have %d and %d entries in HarmonicBondForce, respectively." % (force0.getNumBonds(), force1.getNumBonds())
        n_bonds = force0.getNumBonds()

        dict0, dict1 = {}, {}
        
        i0, i1, r0, k0 = force0.getBondParameters(0)
        unit_r = r0.unit
        unit_k = k0.unit

        for k in range(n_bonds):
            i0, i1, r0, k0 = force0.getBondParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            dict0[i0, i1] = ((r0 / unit_r, k0 / unit_k))

            i0, i1, r0, k0 = force1.getBondParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            dict1[i0, i1] = ((r0 / unit_r, k0 / unit_k))

        assert set(dict0.keys()) == set(dict1.keys()), "Systems have different HarmonicBond Forces"

        for k, parameter_name in enumerate(["r0", "k0"]):
            for (i0, i1) in dict0.keys():
                val0 = dict0[i0, i1][k]
                val1 = dict1[i0, i1][k]
                assert (abs(val0 - val1) / val0) < EPSILON, "Error: Harmonic Bond (%d, %d) has %s values of %f and %f, respectively." % (i0, i1, parameter_name, val0, val1)


    def check_angles(self):
        force0 = self.forces0[1]  # Fix hardcoded lookup later.
        force1 = self.forces1[1]
    
        assert force0.getNumAngles() == force1.getNumAngles(), "Error: Systems have %d and %d entries in HarmonicAngleForce, respectively." % (force0.getNumAngles(), force1.getNumAngles())
        
        n_angles = force0.getNumAngles()

        dict0, dict1 = {}, {}
        
        i0, i1, i2, theta0, k0 = force0.getAngleParameters(0)
        unit_theta = theta0.unit
        unit_k = k0.unit

        for k in range(n_angles):
            i0, i1, i2, theta0, k0 = force0.getAngleParameters(k)
            i0, i1, i2 = reorder_angles(i0, i1, i2)
            dict0[i0, i1, i2] = ((theta0 / unit_theta, k0 / unit_k))

            i0, i1, i2, theta0, k0 = force1.getAngleParameters(k)
            i0, i1, i2 = reorder_angles(i0, i1, i2)
            dict1[i0, i1, i2] = ((theta0 / unit_theta, k0 / unit_k))

        assert set(dict0.keys()) == set(dict1.keys()), "Systems have different HarmonicAngle Forces"

        for k, parameter_name in enumerate(["theta0", "k0"]):
            for (i0, i1, i2) in dict0.keys():
                val0 = dict0[i0, i1, i2][k]
                val1 = dict1[i0, i1, i2][k]
                assert (abs(val0 - val1) / val0) < EPSILON, "Error: Harmonic Angle (%d, %d, %d) has %s values of %f and %f, respectively." % (i0, i1, i2, parameter_name, val0, val1)

    def check_nonbonded(self):
        force0 = self.forces0[0]
        force1 = self.forces1[0]
               
        assert force0.getNumParticles() == force1.getNumParticles(), "Error: Systems have %d and %d particles in NonbondedForce, respectively." % (force0.getNumParticles(), force1.getNumParticles())
        
        n_atoms = force0.getNumParticles()
        
        q, sigma, epsilon = force0.getParticleParameters(0)
        unit_q, unit_sigma, unit_epsilon = q.unit, sigma.unit, epsilon.unit
        
        for k in range(n_atoms):
            q0, sigma0, epsilon0 = force0.getParticleParameters(k)
            q1, sigma1, epsilon1 = force1.getParticleParameters(k)

            q0, sigma0, epsilon0 = q0 / unit_q, sigma0 / unit_sigma, epsilon0 / unit_epsilon
            q1, sigma1, epsilon1 = q1 / unit_q, sigma1 / unit_sigma, epsilon1 / unit_epsilon

            assert (abs(q0 - q1) / q0) < EPSILON, "Error: Particle %d has charges of %f and %f, respectively." % (k, q0, q1)
            if epsilon0 != 0.:
                assert (abs(sigma0 - sigma1) / sigma0) < EPSILON, "Error: Particle %d has sigma of %f and %f, respectively." % (k, sigma0, sigma1)
            else:
                print("Skipping comparison of sigma (%f, %f) on particle %d because epsilon has values %f, %f" % (sigma0, sigma1, k, epsilon0, epsilon1))

            if epsilon0 == 0.:
                denominator = 1.0  # Don't normalize if has value zero
            else:
                denominator = epsilon0  # Normalize to compare relative errors, rather than absolute.

            assert (abs(epsilon0 - epsilon1) / denominator) < EPSILON, "Error: Particle %d has epsilon of %f and %f, respectively." % (k, epsilon0, epsilon1)


        n_exceptions = force0.getNumExceptions()
        assert force0.getNumExceptions() == force1.getNumExceptions(), "Error: Systems have %d and %d exceptions in NonbondedForce, respectively." % (force0.getNumExceptions(), force1.getNumExceptions())


        i0, i1, qq, sigma, epsilon = force0.getExceptionParameters(0)
        unit_qq, unit_sigma, unit_epsilon = qq.unit, sigma.unit, epsilon.unit

        dict0, dict1 = {}, {}
        for k in range(n_exceptions):
            i0, i1, qq, sigma, epsilon = force0.getExceptionParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            dict0[i0, i1] = ((qq / unit_qq, sigma / unit_sigma, epsilon / unit_epsilon))

            i0, i1, qq, sigma, epsilon = force1.getExceptionParameters(k)
            i0, i1 = reorder_bonds(i0, i1)
            dict1[i0, i1] = ((qq / unit_qq, sigma / unit_sigma, epsilon / unit_epsilon))
        
        assert set(dict0.keys()) == set(dict1.keys()), "Systems have different NonBondedForce Exceptions"

        for k, parameter_name in enumerate(["qq", "sigma", "epsilon"]):
            for (i0, i1) in dict0.keys():
                val0 = dict0[i0, i1][k]
                val1 = dict1[i0, i1][k]
                denominator = abs(val0)
                if denominator < EPSILON:
                    denominator = 1.0
                if parameter_name == "sigma" and dict0[i0, i1][2] == 0.0 and dict1[i0, i1][2] == 0.0:
                    continue  # If both epsilon parameters are zero, then sigma doesn't matter so skip the comparison.  
                assert (abs(val0 - val1) / denominator) < EPSILON, "Error: NonBondedForce Exception (%d, %d) has %s values of %f and %f, respectively." % (i0, i1, parameter_name, val0, val1)

    def check_proper_torsions(self):
        bond_force0 = self.forces0[4]
        bond_force1 = self.forces1[4]
        
        force0 = self.forces0[2]
        force1 = self.forces1[2]
        
        bond_set0 = get_symmetrized_bond_set(bond_force0)
        bond_set1 = get_symmetrized_bond_set(bond_force1)
                
        i0, i1, i2, i3, per, phase, k0 = force0.getTorsionParameters(0)
        phase_unit, k0_unit = phase.unit, k0.unit

        dict0, dict1 = {}, {}
        for k in range(force0.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force0.getTorsionParameters(k)
            if not is_proper(i0, i1, i2, i3, bond_set0):
                continue

            i0, i1, i2, i3 = reorder_proper_torsions(i0, i1, i2, i3)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not dict0.has_key((i0, i1, i2, i3)):
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

            if not dict1.has_key((i0, i1, i2, i3)):
                dict1[i0, i1, i2, i3] = []

            dict1[i0, i1, i2, i3].append((per, phase, k0))

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        diff_keys = keys0.symmetric_difference(keys1)

        assert diff_keys == set(), "Systems have different (proper) PeriodicTorsionForce entries: extra keys are: \n%s" % diff_keys

        for (i0, i1, i2, i3) in dict0.keys():
            entries0 = dict0[i0, i1, i2, i3]
            entries1 = dict1[i0, i1, i2, i3]        
            assert len(entries0) == len(entries1), "Error:  (proper) PeriodicTorsionForce entry (%d, %d, %d, %d) has different numbers of terms (%d and %d, respectively)." % (i0, i1, i2, i3, len(entries0), len(entries1))
            
            subdict0 = dict(((per, reduce_precision(phase)), k0) for (per, phase, k0) in entries0)
            subdict1 = dict(((per, reduce_precision(phase)), k0) for (per, phase, k0) in entries1)
            
            assert set(subdict0.keys()) == set(subdict1.keys()), "Error: (proper) PeriodicTorsionForce entry (%d, %d, %d, %d) has different terms." % (i0, i1, i2, i3)
            
            for (per, phase) in subdict0.keys():
                val0 = subdict0[per, phase] 
                val1 = subdict1[per, phase]
                assert (abs(val0 - val1) / val0) < EPSILON, "Error: (proper) PeriodicTorsionForce strength (%d, %d, %d, %d) (%d, %f) has values of %f and %f, respectively." % (i0, i1, i2, i3, per, phase, val0, val1)

    def check_improper_torsions(self):
        bond_force0 = self.forces0[4]
        bond_force1 = self.forces1[4]
        
        force0 = self.forces0[2]
        force1 = self.forces1[2]
        
        bond_set0 = get_symmetrized_bond_set(bond_force0)
        bond_set1 = get_symmetrized_bond_set(bond_force1)
                
        i0, i1, i2, i3, per, phase, k0 = force0.getTorsionParameters(0)
        phase_unit, k0_unit = phase.unit, k0.unit

        dict0, dict1 = {}, {}
        for k in range(force0.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force0.getTorsionParameters(k)
            if is_proper(i0, i1, i2, i3, bond_set0):
                continue

            i0, i1, i2, i3 = reorder_improper_torsions(i0, i1, i2, i3)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not dict0.has_key((i0, i1, i2, i3)):
                dict0[i0, i1, i2, i3] = []

            dict0[i0, i1, i2, i3].append((per, phase, k0))

        for k in range(force1.getNumTorsions()):
            i0, i1, i2, i3, per, phase, k0 = force1.getTorsionParameters(k)
            if is_proper(i0, i1, i2, i3, bond_set1):
                continue

            i0, i1, i2, i3 = reorder_improper_torsions(i0, i1, i2, i3)

            phase, k0 = phase / phase_unit, k0 / k0_unit
            if k0 == 0.0:
                continue

            if not dict1.has_key((i0, i1, i2, i3)):
                dict1[i0, i1, i2, i3] = []

            dict1[i0, i1, i2, i3].append((per, phase, k0))

        self.dict0, self.dict1 = dict0, dict1

        keys0 = set(dict0.keys())
        keys1 = set(dict1.keys())
        diff_keys = keys0.symmetric_difference(keys1)

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
                assert (abs(val0 - val1) / val0) < EPSILON, "Error: (improper) PeriodicTorsionForce strength (%d, %d, %d, %d) (%d, %f) has values of %f and %f, respectively." % (i0, i1, i2, i3, per, phase, val0, val1)

