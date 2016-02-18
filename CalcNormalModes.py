#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Mòdul per calcular NormalModes a partir d'una molècula de chimera utilitzant ProDy
"""

import chimera
import numpy
import prody


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


def alg1(moldy, n=7):
    """
    Coarse Grain Algorithm 1: groups per residues
    """
    group = 1
    for chain in moldy.iterChains():
        selection = moldy.select('chain {}'.format(chain.getChid()))
        num_residues = selection.getResnums()[-1]
        #num_residues = selection.getHierView().numResidues()
        for start, end in chunker(num_residues, n):
            try:
                moldy.select('chain {} and resnum {} to {}'.format(
                    chain.getChid(), start, end)).setBetas(group)
                group += 1
            except:  # afegir tipus d'error
                pass
    return moldy


def alg2(moldy, n=100):
    """
    Coarse Grain Algorithm 2: groups per mass percentage
    n: tant per cent of total mass per group
    """

    group = 1

    M = sum(moldy.getMasses())
    m = M/n
    mass = None

    for chain in moldy.iterChains():
        selection = moldy.select('chain {}'.format(chain.getChid()))
        num_residues = selection.getResnums()[-1]
        num_atoms = selection.numAtoms()
        mass = 0.

        for atom in iter(selection):  # atoms a la cadena? residus?
            atom.setBeta(group)
            mass += atom.getMass()
            if mass > m:
                mass = 0.
                group += 1
        group += 1
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
    return calc_normal_modes(mol, alg1)
