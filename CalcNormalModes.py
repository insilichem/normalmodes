#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calc Normal Modes from a chimera molecule using prody
"""

import chimera
import numpy
import prody
import CoarseGrainAlgorithm as CGAlg
from StructBio.EnergyFunctions import LJ

PRODY2CHIMERA = {}
ENERGY = None
MOLECULE = None
EPSILON = 0.15
SIGMA = 4.0

SIGMA_6 = SIGMA**6
SIGMA_12 = SIGMA_6**2


def calc_normal_modes(mol, coarse_grain=None,
                      n_algorithm=None, n_modes=None,
                      LJ=False, mass_weighted=True):
    """
    Parameters
    ----------
    mol : chimera.Molecule or prody.AtomGroup
    CoarseGrain : callable, optional, default=None
        coarseGrain(prm) wich make mol.select().setBetas(i) where i
        is the index Coarse Grain group
        Where prm is prody AtomGroup

    Returns
    -------
        modes ProDy like ANM or RTB
    """
    global ENERGY, MOLECULE
    MOLECULE = mol
    if isinstance(mol, chimera.Molecule):
        delete_waters(mol)
        moldy = chimera2prody(mol)
    elif isinstance(mol, prody.AtomGroup):
        moldy = mol
    else:
        pass  # Escriure error

    modes = None
    if coarse_grain == ('Residues' or 'Mas'):
        if n_algorithm > 0:
            moldy = coarse_grain(moldy, n_algorithm)
        else:
            moldy = coarse_grain(moldy)
        title = 'normal modes for {} using algorithm {}'.format(
            moldy.getTitle(), coarse_grain.title)
        modes = prody.RTB(title)
        if LJ:
            modes.buildHessian(moldy, blocks=moldy.getBetas(), gamma=gamma_lennardjones)
        else:
            modes.buildHessian(moldy.getCoords(), blocks=moldy.getBetas())
        
    else:
        modes = prody.ANM('normal modes for {}'.format(moldy.getTitle()))
        if LJ:
            modes.buildHessian(moldy, gamma=gamma_lennardjones)
        else:
            modes.buildHessian(moldy)

    if mass_weighted:
        N = moldy.numAtoms()
        masses = numpy.zeros(3*N,float)
        for i in xrange(len(moldy.getMasses())):
            masses[i*3:i*3+3] = moldy.getMasses()[i]
        reduced_masses = 1/numpy.sqrt(masses)
        v, H = reduced_masses, modes.getHessian()
        hessian = ((v*H).T*v).T
        modes.setHessian(hessian)

    modes.calcModes(n_modes=n_modes)
    return modes


def chimera2prody(mol):
    """
    Function that transforms a chimera molecule into a prody atom group

    Parameters
    ----------
    mol: chimera.Molecule

    Returns
    -------
    moldy: prody.AtomGroup()
    """
    moldy = prody.AtomGroup()
    try:
        # moldy.setCoords(numpy.array([tuple(atm.coord()) for atm in mol.atoms], dtype=float))
        # moldy.setElements([atm.element.name for atm in mol.atoms])
        # moldy.setNames([atm.name for atm in mol.atoms])
        # moldy.setResnums([atm.residue.id.position for atm in mol.atoms])
        # moldy.setChids([atm.residue.id.chainId for atm in mol.atoms])
        # moldy.setBetas([atm.bfactor for atm in mol.atoms])
        # moldy.setMasses([atm.element.mass for atm in mol.atoms])
        # # moldy.setBonds([[bond.atoms[0].coordIndex, bond.atoms[1].coordIndex]for bond in mol.bonds])
        # moldy.setTitle(str(mol.name))

        # n = len(mol.atoms)
        # coords = numpy.zeros(n)
        coords, elements, names, resnums, chids, betas, masses = \
            [], [], [], [], [], [], []
        e = {}
        global PRODY2CHIMERA
        offset_chimera_residue = min(r.id.position for r in mol.residues)
        for i, atm in enumerate(mol.atoms):
            PRODY2CHIMERA[i] = atm.coordIndex
            e[atm.coordIndex] = i
            coords.append(tuple(atm.coord()))  # array documentation to improve
            elements.append(atm.element.name)
            names.append(atm.name)
            resnums.append(atm.residue.id.position - offset_chimera_residue)
            chids.append(atm.residue.id.chainId)
            masses.append(atm.element.mass)
            betas.append(atm.bfactor)

        moldy.setCoords(coords)
        moldy.setElements(elements)
        moldy.setNames(names)
        moldy.setResnums(resnums)
        moldy.setChids(chids)
        moldy.setBetas(betas)
        moldy.setMasses(masses)

        moldy.setBonds([[e[bond.atoms[0].coordIndex], e[bond.atoms[1].coordIndex]]
                        for bond in mol.bonds])

        moldy.setTitle(mol.name)

    except AttributeError:
        raise AttributeError('mol must be a chimera.Molecule')
    return moldy


# def mass(*items):
#     try:
#         return sum([items])
#     if isinstance(items[0], chimera.Atom):
#         return sum([a.element.mass for a in items])
#     elif isinstance(items[0], chimera.Residue):
#         return sum([a.element.mass for r in items for a in r.atoms])

    # try:
    #     return sum([a.element.mass for a in items])
    # except AttributeError:
    #     return sum([a.element.mass for r in items for a in r.atoms])


def translation(vector, point):
    t1 = [vector[0]/point[0], point[0]/(2*point[1]), point[0]/(2*point[2])]
    t2 = [point[1]/(2*point[0]), vector[1]/point[1], point[1]/(2*point[2])]
    t3 = [point[2]/(2*point[0]), point[2]/(2*point[1]), vector[2]/point[2]]
    translation = numpy.array((t1, t2, t3))
    return translation


def lennard_jones(molecule):
    """
    molecule must be a chimera molecule
    return:
        a matrix between pairs of chimera indices and its lennard Jones Potential
    """
    LJ_tools = LJ.LJ_evaluation()
    n_atoms = len(molecule.atoms)
    energy = numpy.zeros((n_atoms,n_atoms),float)
    for i in xrange(n_atoms):
        for j in xrange(n_atoms):
            energy[i,j] = LJ_tools.atomPairEnergy(molecule.atoms[i],molecule.atoms[j])
            if energy[i,j] < 1.0:
                energy[i,j] = 1.0
    
    return energy

def gamma_lennardjones(dist2, i, j):
    """
    given a pair of atom indices i and j, return the lennard jones energy
    """
    # energy = 1.0
    i, j = PRODY2CHIMERA[i], PRODY2CHIMERA[j]
    atom_i, atom_j = MOLECULE.atoms[i], MOLECULE.atoms[j]
    r = numpy.linalg.norm(numpy.array(atom_i.coord()) - numpy.array(atom_j.coord()))
    r_6 = r**6
    r_12 = r_6**2

    force = -4.0*EPSILON*(SIGMA_12/r_12-SIGMA_6/r_6)
    force *= dist2

    # if force > energy:
    return min(force,100)
    # else:
    #     return energy


def sample_and_translation(moldy, modes):
    """
    calculatea new distribution
    make the hessian matrix: with prody and with linear-transformation
    """
    new_molecule = prody.AtomGroup()
    ensemble = prody.sampleModes(modes=modes, atoms=moldy, n_confs=1)
    new_molecule.setCoords(ensemble.getCoords())
    vector = new_molecule.getCoords()-moldy.getCoords()
    num_atoms = new_molecule.numAtoms()
    T = numpy.zeros((3*num_atoms, 3*num_atoms))
    for i in xrange(num_atoms):
        point = moldy.getCoords()[i]
        vector = new_molecule.getCoords()[i]-point
        T[i*3:(i+1)*3, i*3:(i+1)*3] = translation(vector, point)
    new_hessian = numpy.transpose(T)*modes.getCoordsHessian()*T

    new_modes = calc_normal_modes(new_molecule)

    return new_modes, new_hessian


def main(mol):
    return calc_normal_modes(mol, CGAlg.alg1)

def delete_waters(molecule):
    for atom in molecule.atoms:
        if atom.residue.type == 'HOH':
            molecule.deleteAtom(atom)