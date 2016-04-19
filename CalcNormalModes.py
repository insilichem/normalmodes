#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calc Normal Modes from a chimera molecule using prody
"""

import chimera
import numpy
import prody
import CoarseGrainAlgorithm as CGAlg


def calc_normal_modes(mol, coarse_grain=None, n_algorithm=None, n_modes=None):
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
    if isinstance(mol, chimera.Molecule):
        moldy = chimera2prody(mol)
    elif isinstance(mol, prody.AtomGroup):
        moldy = mol
    else:
        pass  # Escriure error

    modes = None
    if coarse_grain:
        if n_algorithm > 0:
            moldy = coarse_grain(moldy, n_algorithm)
        else:
            moldy = coarse_grain(moldy)
        title = 'normal modes for {} using algorithm {}'.format(
            moldy.getTitle(), coarse_grain.title)
        modes = prody.RTB(title)
        modes.buildHessian(moldy.getCoords(), moldy.getBetas())
        modes.calcModes(n_modes=n_modes)
    else:
        modes = prody.ANM('normal modes for {}'.format(moldy.getTitle()))
        modes.buildHessian(moldy)
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
        d, e = {}, {}
        offset_chimera_residue = min(r.id.position for r in mol.residues)
        for i, atm in enumerate(mol.atoms):
            d[i] = atm.coordIndex
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
