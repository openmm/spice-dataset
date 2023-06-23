import openmm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit.topology import Molecule
from openff.units import unit as ffunit
import numpy as np
import h5py

# Create a box of water and model it with AMOEBA.

modeller = app.Modeller(app.Topology(), [])
modeller.addSolvent(app.ForceField('tip3pfb.xml'), boxSize=(2.0, 2.0, 2.0)*unit.nanometers)
forcefield = app.ForceField('amoeba2018.xml')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.7*unit.nanometers, vdwCutoff=0.9*unit.nanometers, mutualInducedTargetEpsilon=1e-5, polarization='mutual')
integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()

# Simulate it for 10 ns, saving a state every 10 ps.

states = []
for i in range(1000):
    simulation.step(10000)
    states.append(simulation.context.getState(getPositions=True, enforcePeriodicBox=True))

# For each state, find the 30 water molecules closest to the center of the box.

oxygen = [atom.index for atom in simulation.topology.atoms() if atom.element == app.element.oxygen]
clusters = []
for state in states:
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer) - np.array([1.0, 1.0, 1.0])
    dist2 = np.sum(pos[oxygen]*pos[oxygen], axis=1)
    closest = np.argsort(dist2)[:30]
    indices = []
    for i in closest:
        indices += [oxygen[i], oxygen[i]+1, oxygen[i]+2]
    clusters.append(pos[indices])

# Use OpenFF to canonically order the atoms and generate the SMILES string.

mol = Molecule()
for i in range(30):
    mol.add_atom(app.element.oxygen.atomic_number, 0, False)
    mol.add_atom(app.element.hydrogen.atomic_number, 0, False)
    mol.add_atom(app.element.hydrogen.atomic_number, 0, False)
    mol.add_bond(3*i, 3*i+1, 1, False)
    mol.add_bond(3*i, 3*i+2, 1, False)
mol._conformers = None
for c in clusters:
    mol.add_conformer(c*unit.nanometers)
mol = mol.canonical_order_atoms()
smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
conformations = [c.m_as(ffunit.nanometers) for c in mol.conformers]

# Save the clusters to a HDF5 file.

outputfile = h5py.File('water.hdf5', 'w')
group = outputfile.create_group('water')
group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
ds.attrs['units'] = 'nanometers'
