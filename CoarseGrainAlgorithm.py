#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Algorithm to do Coarse Grain
Same beta for same group
To pass to CalcNormalModes.calc_normal_modes()
"""

import prody


def alg1(moldy, n=7):
    """
    Coarse Grain Algorithm 1: groups per residues

    Parameters
    ----------
    moldy : prody.AtomGroup
    n : int, optional, default=7
        number of residues per group

    Returns
    ------
    moldy: prody.AtomGroup
        New betas added
    """
    group = 1
    for chain in moldy.iterChains():
        selection = moldy.select('chain {}'.format(chain.getChid()))
        num_residues = selection.getResnums()
        for a, b in chunker(len(num_residues), n):
            try:
                chain, start, end = chain.getChid(), num_residues[a], num_residues[b]
                print chain, a, ':', start, b, ':', end
                selector = 'chain {} and resnum {} to {}'.format(chain, start, end)
                selection = moldy.select(selector)
                selection.setBetas(group)
                group += 1
            except AttributeError as e:
                print('Warning: {}'.format(e))
                pass
    return moldy


def alg2(moldy, n=100):
    """
    Coarse Grain Algorithm 2: groups per mass percentage

    Parameters
    ----------
    moldy : prody.AtomGroup
    n : int, optional, default=100
        number of groups

    Returns
    -------
    moldy: prody.AtomGroup
        New Betas added
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
