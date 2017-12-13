#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from threading import Thread
from Queue import Queue
# import copy
import chimera
from chimera.tasks import Task
import prody
from cclib.parser import Gaussian
from NormalModesTable import NormalModesTableDialog
# from new_gui import NormalModesResultsDialog, NormalModesMovieDialog

ALGORITHMS = {'Full atom': 'full_atom',
              'Group by residues': 'residues',
              'Group by mass': 'mass',
              'Group by graph': 'graph',
              'Extend from C-alpha': 'calpha'}


class Controller(object):

    """docstring for Controller"""


    def __init__(self, gui=None, *args, **kwargs):
        self.gui = gui
        self.vibrations = None
        self._molecules = None
        self.ENGINES = dict(prody=self._run_prody,
                            gaussian=self._run_gaussian)

    def run(self):
        self._failure()
        self.gui.buttonWidgets['Run']['state'] = 'disabled'
        self.gui.buttonWidgets['Run']['text'] = 'Running...'
        engine = self.ENGINES[self.gui.engine]
        # Run external task with the following to prevent UI freezes
        task = StatusTask('Normal Modes Analysis', cancelCB=self._failure)
        task.put = task.updateStatus
        queue = ChimeraTaskQueue(task=task)
        engine(queue=queue)
        task.put('Done!')
        self.gui.buttonWidgets['Run']['state'] = 'normal'
        self.gui.buttonWidgets['Run']['text'] = 'Run'
        self._success()
        task.finished()

    def _run_prody(self, queue=None, threaded=True):
        self._molecules = self.gui.ui_molecules.getvalue()
        if not self._molecules:
            self.gui.buttonWidgets['Run']['state'] = 'normal'
            self.gui.buttonWidgets['Run']['text'] = 'Run'
            raise chimera.UserError("Please select at least one molecule")

        algorithm = ALGORITHMS[self.gui.ui_algorithms_menu.getvalue()]
        algorithm_param = int(self.gui.ui_algorithms_param.get())
        n_modes = int(self.gui.ui_n_modes.get())
        cutoff = float(self.gui.ui_cutoff.get())
        gamma = float(self.gui.ui_gamma.get())
        mw = 'Mass-Weighted' in self.gui.ui_extra_options_chk.getvalue()
        if threaded:
            thread = Thread(target=VibrationalMolecule.from_chimera,
                            args=(self._molecules,),
                            kwargs=dict(algorithm=algorithm,
                                        n=algorithm_param,
                                        queue=queue,
                                        max_modes=n_modes,
                                        gamma=gamma,
                                        mass_weighted=mw,
                                        cutoff=cutoff))
            thread.start()
            while thread.isAlive():
                chimera.tkgui.app.update()
            self.vibrations = queue.get_nowait()
        else:
            self.vibrations = VibrationalMolecule.from_chimera(self._molecules,
                                                              algorithm=algorithm,
                                                              n=algorithm_param,
                                                              queue=queue.task,
                                                              max_modes=n_modes,
                                                              gamma=gamma,
                                                              mass_weighted=mw,
                                                              cutoff=cutoff)
        return True

    def _run_gaussian(self, **kwargs):
        path = self.gui.ui_gaussian_file_entry.get()
        if not os.path.exists(path):
            raise chimera.UserError("File {} not available".format(path))
        self.vibrations = VibrationalMolecule.from_gaussian(path)
        self._molecules = self.vibrations.molecule
        return True

    def _success(self):
        if self.vibrations is None:
            return
        frequencies = self.vibrations.frequencies()[1]
        # results_dialog = NormalModesResultsDialog(self.gui, controller=self)
        # results_dialog.enter()
        # results_dialog.populate_data(frequencies)
        # movie_dialog = NormalModesMovieDialog(self.gui, controller=self)
        # movie_dialog.enter()
        dialog = NormalModesTableDialog(self._molecules, self.vibrations.modes)

    def _failure(self):
        self.vibrations = None
        # print some kind of message


class VibrationalMolecule(object):
    """
    class to handle prody normal modes
    """

    def __init__(self, molecule, modes, task=None):
        self.molecule = molecule
        self.modes = modes

    @classmethod
    def from_chimera(cls, chimera_molecule, queue=None, **kwargs):
        """
        initializes the VibrationalMolecule class with a chimera molecule
        """
        molecule, chimera2prody = convert_chimera_molecule_to_prody(chimera_molecule)
        modes = calculate_vibrations(molecule, queue=queue, **kwargs)
        instance = cls(chimera_molecule, modes)
        if queue is not None:
            queue.put_nowait(instance)
        return instance

    @classmethod
    def from_gaussian(cls, gaussian_path, task=None):
        """
        initializes from a gaussian output file
        """
        gaussian_parser = Gaussian(gaussian_path).parse()
        modes = gaussian_modes(gaussian_parser)
        chimera_molecule = convert_gaussian_molecule_to_chimera(gaussian_parser)
        chimera.openModels.add([chimera_molecule])
        return cls(chimera_molecule, modes)

    def frequencies(self):
        """
        normal modes vectors and frequencies
        """
        return self.modes.getEigvecs(), self.modes.getEigvals()

    def trajectory(self, mode, n_steps, rmsd=1.5):
        """
        given a number of a selected mode and a number of stemps returns
        a list with n_steps coordsets
        """
        ensemble = prody.traverseMode(mode=mode, atoms=self.molecule,
                                      n_steps=n_steps, rmsd=rmsd)
        return ensemble.getCoordsets()

    def sample(self,modes,rmsd):
        """
        modes : list
        rmsd : float
        """
        conformation = prody.sampleModes(modes=modes, atoms=self.molecule,
                                         n_confs=1, rmsd=rmsd)
        return conformation.getCoords()


class ChimeraTaskQueue(Queue):

    def __init__(self, task=None, *args, **kwargs):
        self.task = task
        Queue.__init__(self, *args, **kwargs)

    def _put(self, item):
        if isinstance(item, basestring) and self.task is not None:
            self.task.updateStatus(item)
        else:
            Queue._put(self, item)


class StatusTask(Task):

    def updateStatus(self, msg):
        chimera.statusline.show_message('{}: {}'.format(self.title, msg),
            clickCallback=self._status_click_callback)
        Task.updateStatus(self, msg)

    def _status_click_callback(self, *a, **kw):
        chimera.dialogs.display(chimera.tasks.TaskPanel.name)


def convert_chimera_molecule_to_prody(molecule):
    """
    Function that transforms a chimera molecule into a prody atom group

    Parameters
    ----------
    molecule : chimera.Molecule

    Returns
    -------
    prody_molecule : prody.AtomGroup()
    chimera2prody : dict
        dictionary: chimera2prody[chimera_atom.coordIndex] = i-thm element prody getCoords() array
    """
    if not isinstance(molecule, chimera.Molecule):
        raise TypeError('`molecule` must be a chimera.Molecule')

    prody_molecule = prody.AtomGroup()
    coords, elements, serials = [], [], []
    names, resnums, resnames = [], [], []
    chids, betas, masses = [], [], []
    chimera2prody = {}
    offset_chimera_residue = min(r.id.position for r in molecule.residues)

    for i, atm in enumerate(molecule.atoms):
        chimera2prody[atm.serialNumber] = i
        coords.append(tuple(atm.coord()))  # array documentation to improve
        elements.append(atm.element.name)
        serials.append(atm.serialNumber)
        names.append(atm.name)
        resnums.append(atm.residue.id.position - offset_chimera_residue)
        resnames.append(atm.residue.type)
        chids.append(atm.residue.id.chainId)
        masses.append(atm.element.mass)
        betas.append(atm.bfactor)

    prody_molecule.setCoords(coords)
    prody_molecule.setElements(elements)
    prody_molecule.setSerials(serials)
    prody_molecule.setNames(names)
    prody_molecule.setResnums(resnums)
    prody_molecule.setResnames(resnames)
    prody_molecule.setChids(chids)
    prody_molecule.setBetas(betas)
    prody_molecule.setMasses(masses)
    prody_molecule.setTitle(str(molecule.name))
    prody_molecule.setBonds([(chimera2prody[bond.atoms[0].serialNumber],
                                chimera2prody[bond.atoms[1].serialNumber]) for bond in molecule.bonds])

    return prody_molecule, chimera2prody


def convert_gaussian_molecule_to_chimera(path_or_parser):
    if isinstance(path_or_parser, basestring):
        parsed = Gaussian(path_or_parser).parse()
    else:
        parsed = path_or_parser

    elements = parsed.atomnos
    coords = parsed.atomcoords[-1]
    m = chimera.Molecule()
    r = m.newResidue("UNK", " ", 1, " ")
    for i, (element, xyz) in enumerate(zip(elements, coords)):
        element = chimera.Element(element)
        a = m.newAtom('{}{}'.format(element, i), element)
        r.addAtom(a)
        a.setCoord(chimera.Point(*xyz))
        a.serialNumber = i
    chimera.connectMolecule(m)
    return m


def calculate_vibrations(molecule, max_modes=20, algorithm='calpha', queue=None,
                         mass_weighted=False, cutoff=15.0, gamma=1.0,
                         **options):
    """
    Parameters
    ----------
    molecule : prody.AtomGroup
    nax_modes : int
        number of modes to calculate
    algorithm : callable, optional, default=None
        coarseGrain(prm) which make molecule.select().setBetas(i) where i
        is the index Coarse Grain group
        and prm is prody AtomGroup
    options : dict, optional
        Parameters for algorithm callable

    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    if queue is None:
        queue = Queue()
    modes = None
    hessian_kwargs = dict(cutoff=cutoff, gamma=gamma)
    if algorithm in ('residues', 'mass', 'graph'):
        queue.put('Building model...')
        title = 'normal modes for {}'.format(molecule.getTitle())
        molecule = GROUPERS[algorithm](molecule, **options)
        modes = prody.RTB(title)
        queue.put('Building hessian...')
        modes.buildHessian(molecule.getCoords(), molecule.getBetas(), **hessian_kwargs)
        if mass_weighted:
            queue.put('Mass weighting...')
            _mass_weighted_hessian(molecule, modes)
        queue.put('Calculating {} modes...'.format(max_modes))
        modes.calcModes(n_modes=max_modes)
    elif algorithm == 'calpha':
        queue.put('Building model...')
        calphas_modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        calphas = molecule = molecule.select(algorithm)
        queue.put('Building hessian...')
        calphas_modes.buildHessian(calphas, **hessian_kwargs)
        if mass_weighted:
            queue.put('Mass weighting...')
            _mass_weighted_hessian(molecule, modes)
        queue.put('Calculating {} modes...'.format(max_modes))
        calphas_modes.calcModes(n_modes=max_modes)
        queue.put('Extending model...')
        modes = prody.extendModel(calphas_modes, calphas, molecule, norm=True)[0]
    else:
        queue.put('Building model...')
        modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        queue.put('Building hessian...')
        modes.buildHessian(molecule, **hessian_kwargs)
        if mass_weighted:
            queue.put('Mass weighting...')
            _mass_weighted_hessian(molecule, modes)
        queue.put('Calculating {} modes...'.format(max_modes))
        modes.calcModes(n_modes=max_modes)
    print(modes)
    return modes


def _mass_weighted_hessian(molecule, modes):
    N = molecule.numAtoms()
    masses = np.zeros(3 * N, float)
    for i in xrange(len(molecule.getMasses())):
        masses[i * 3:i * 3 + 3] = molecule.getMasses()[i]
    reduced_masses = 1 / np.sqrt(masses)
    v, H = reduced_masses, modes.getHessian()
    hessian = ((v * H).T * v).T
    modes.setHessian(hessian)


def gaussian_modes(path_or_parser):
    """
    Read the modes
    Create a prody.modes instance

    Parameters
    ----------
    path_or_parser : str or cclib.Gaussian
        gaussian frequencies output path
        or already parsed cclib.Gaussian instance

    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    if isinstance(path_or_parser, basestring):
        parsed = Gaussian(path_or_parser).parse()
    else:
        parsed = path_or_parser
    shape = parsed.vibdisps.shape
    modes_vectors = parsed.vibdisps.reshape(shape[0], shape[1]*shape[2]).T
    modes_frequencies = np.abs(parsed.vibfreqs)
    modes = prody.NMA()
    modes.setEigens(vectors=modes_vectors, values=modes_frequencies)
    return modes


def group_by_residues(molecule, n=7):
    """
    Coarse Grain Algorithm 1: groups per residues

    Parameters
    ----------
    molecule : prody.AtomGroup
    n : int, optional, default=7
        number of residues per group

    Returns
    ------
    molecule : prody.AtomGroup
        New betas added
    """
    group = 1
    for chain in molecule.iterChains():
        residues_indices = sorted(list(set(chain.getResnums())))
        chain_name = chain.getChid()
        for a, b in chunker(len(residues_indices), n):
            try:
                start, end = residues_indices[a-1], residues_indices[b-1]
                selector = 'chain {} and resnum {} to {}'.format(chain_name, start, end)
                selection = molecule.select(selector)
                selection.setBetas(group)
                group += 1
            except AttributeError as e:
                print(e)
    return molecule


def group_by_mass(molecule, n=100):
    """
    Coarse Grain Algorithm 2: groups per mass percentage

    Parameters
    ----------
    molecule : prody.AtomGroup
    n: int, optional, default=100
        Intended number of groups. The mass of the system will be divided by this number,
        and each group will have the corresponding proportional mass. However, the final
        number of groups can be slightly different.

    Returns
    -------
    molecule: prody.AtomGroup
        New Betas added
    """
    group = 1

    total_mass = sum(molecule.getMasses())
    chunk_mass = total_mass/n

    for chain in molecule.iterChains():
        selection = molecule.select('chain {}'.format(chain.getChid()))
        mass_accumulator = 0.

        for atom in iter(selection):
            atom.setBeta(group)
            mass_accumulator += atom.getMass()
            if mass_accumulator > chunk_mass:
                mass_accumulator = 0.
                group += 1
        group += 1
    return molecule


def group_by_graph(moldy, n=2):
    """
    Coarse Grain Algorithm 3: Graph algorithm.
        New group when a vertice: have more than n,
                                  have 0 edges
                                  new chain

    Parameters
    ----------
    moldy : prody.AtomGroup
    n : int, optional, default=2
        maximum bonds number

    Returns
    -------
    moldy: prody.AtomGroup
        New Betas added
    """
    group = 1
    title = 'Graphs'

    for chain in moldy.iterChains():
        # selection = moldy.select('chain {}'.format(chain.getChid()))
        # for atom in iter(selection):
        for atom in chain.iterAtoms():
            atom.setBeta(group)
            if atom.numBonds() > n or atom.numBonds() == 0:
                group += 1
        group += 1
    return moldy


def chunker(end, n):
    """
    divide end integers in closed groups of n
    """
    for i in range(0, end-n+1, n):
        yield i+1, i+n
    if end % n:
        yield end-end % n+1, end


GROUPERS = dict(residues= group_by_residues,
               mass=group_by_mass,
               graph=group_by_graph)
