#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calc Normal Modes from a chimera molecule using prody
"""

import chimera
import numpy
import prody
import coarseGrainAlgorithm as CGAlg


def calc_normal_modes(mol, coarse_grain=None,):
    """
    Parameters
    ----------
    mol : chimera.Molecule
    CoarseGrain : callable
        coarseGrain(prm) wich make mol.select().setBetas(i) where i
        is the index Coarse Grain group
        Where prm is prody AtomGroup

    Returns
    -------
        modes ProDy like ANM or RTB
    """

    moldy = chimera2prody(mol)
    # Descomentar la següent línia si la molecula mol ja és de prody
    # moldy = mol
    modes = None
    if coarse_grain:
        moldy = coarse_grain(moldy)
        modes = prody.RTB('RTB for {} using algorithm {}'.format(moldy.getTitle, 1))
        modes.buildHessian(moldy.getCoords(), moldy.getBetas())
        modes.calcModes()
    else:
        modes = prody.ANM('ANM for {}'.format(moldy.getTitle))
        modes.buildHessian(moldy)
        modes.calcModes()
    return modes


def chimera2prody(mol):
    moldy = prody.AtomGroup()
    try:
        moldy.setCoords(numpy.array([tuple(atm.coord()) for atm in mol.atoms], dtype=float))
        moldy.setElements([atm.element.name for atm in mol.atoms])
        moldy.setNames([atm.name for atm in mol.atoms])
        moldy.setResnums([atm.residue.id.position for atm in mol.atoms])
        moldy.setChids([atm.residue.id.chainId for atm in mol.atoms])
        moldy.setBetas([atm.bfactor for atm in mol.atoms])
        moldy.setMasses([atm.element.mass for atm in mol.atoms])
        moldy.setTitle(str(mol.name))
        """
        n = len(mol.atoms)
        coords, elements, names, resnums, masses = numpy.zeros(n), [], [], [], []
        for atm in mol.atoms:
            coords.append(tuple(atm.coord()))#array documentation
            elements.append(atm.element.name)
            names.append(atm.name)
            resnums.append(atm.residue.id.position)
            masses.append(atm.element.mass)
        moldy.setCoords(coords)
        """

    except AttributeError:
        raise AttributeError('mol must be a chimera.Molecule')
    return moldy


def chunker(end, n):
    for i in range(0, end-n+1, n):
        yield i+1, i+n
    if end % n:
        yield end-end % n+1, end


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


def main(mol):
    return calc_normal_modes(mol, CGAlg.alg1)
