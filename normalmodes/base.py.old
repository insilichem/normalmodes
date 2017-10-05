#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from calc_normal_modes import calc_normal_modes
import algorithms
from NormalModesTable import NormalModesTableDialog


class nmod(object):

    ALGORITHMS = dict(Residues=algorithms.alg1,
                      Mas=algorithms.alg2)
                      
    def __init__(self, mols, proc, prodyalgorithm=None, n_algorithm=None,
                 LJ=False, mass_weighted=True):
        self.molecules = mols
        algorithm = self.ALGORITHMS.get(prodyalgorithm)
        self.modes = calc_normal_modes(self.molecules[0], algorithm, n_algorithm,
                                       LJ=LJ, mass_weighted=mass_weighted)
        self.visualization()
        self.modes = nmod

    def visualization(self):
        self.freq_dialog = NormalModesTableDialog(self.molecules[0], self.modes)


class gaussian_nmod(object):

    def __init__(self, mols, proc, path):
        self.molecules = mols
        self.modes = gaussian2prody(path)
        self.visualization()

    def visualization(self):
        NormalModesTableDialog(self.molecules[0], self.modes)


def gaussian2prody(path):
    return
