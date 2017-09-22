
import os
import numpy as np
import prody
# import copy
import chimera
from cclib.parser import Gaussian
from NormalModesTable import NormalModesTableDialog
# from new_gui import NormalModesResultsDialog, NormalModesMovieDialog


class Controller(object):

    """docstring for Controller"""


    def __init__(self, gui=None, *args, **kwargs):
        self.gui = gui
        self.vibrations = None
        self._molecules = None
        self.ENGINES = dict(prody=self._run_prody,
                            gaussian=self._run_gaussian)

    def run(self):
        engine = self.ENGINES[self.gui.engine]
        # Run external task with the following to prevent UI freezes
        ok = engine()
        if ok:
            self._success()
        else:
            self._failure()
        

    def _run_prody(self):
        self._molecules = self.gui.ui_molecules.getvalue()
        if not self._molecules:
            raise chimera.UserError("Please select at least one molecule")
        
        algorithm = self.gui.ui_algorithms_menu.getvalue()
        algorithm_param = self.gui.ui_algorithms_param.get()
        extra_option = self.gui.ui_extra_options_chk.getvalue()
        self.vibrations = VibrationalMolecule.from_chimera(self._molecules, 
                                                           algorithm=algorithm, 
                                                           n=algorithm_param)
        return True

    def _run_gaussian(self):
        path = self.gui.ui_gaussian_file_entry.get()
        if not os.path.exists(path):
            raise chimera.UserError("File {} not available".format(path))
        self.vibrations = VibrationalMolecule.from_gaussian(path)
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

    def __init__(self, molecule, modes):
        self.molecule = molecule
        self.modes = modes

    @classmethod
    def from_chimera(cls, chimera_molecule, **kwargs):
        """
        initializes the VibrationalMolecule class with a chimera molecule
        """
        molecule, chimera2prody = convert_chimera_molecule_to_prody(chimera_molecule)
        modes = calculate_vibrations(molecule, **kwargs)
        return cls(molecule, modes)

    @classmethod
    def from_gaussian(cls, gaussian_path):
        """
        initializes from a gaussian output file
        """
        gaussian_parser = Gaussian(gaussian_path).parse()
        molecule = prody.AtomGroup()
        molecule.setCoords(gaussian_parser.atomcoords)
        shape = gaussian_parser.vibdisps.shape
        modes_vectors = gaussian_parser.vibdisps.reshape(shape[0], shape[1]*shape[2]).T
        modes_frequencies = np.abs(gaussian_parser.vibfreqs)
        modes = prody.NMA()
        modes.setEigens(vectors=modes_vectors, values=modes_frequencies)
        return cls(molecule, modes)

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
        ensemble = prody.traverseMode(mode=mode,atoms=self.molecule,
                                      n_steps=n_steps,rmsd=rmsd)
        return ensemble.getCoordsets()

    def sample(self,modes,rmsd):
        """
        modes : list
        rmsd : float
        """
        conformation = prody.sampleModes(modes=modes, atoms=self.molecule,
                                         n_confs=1, rmsd=rmsd)
        return conformation.getCoords()


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
    prody_molecule = prody.AtomGroup()
    try:
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

    except AttributeError:
        raise TypeError('Attribute not found. Molecule must be a chimera.Molecule')

    return prody_molecule, chimera2prody


def calculate_vibrations(molecule, max_modes=20, algorithm='calpha', **options):
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
    modes = None
    if algorithm in ['residues', 'mass']:
        title = 'normal modes for {}'.format(molecule.getTitle())
        molecule = algorithm(molecule, **options)
        modes = prody.RTB(title)
        modes.buildHessian(molecule.getCoords(), molecule.getBetas())
        modes.calcModes(n_modes=max_modes)
    elif algorithm == 'calpha':
        calphas_modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        calphas = molecule = molecule.select(algorithm)
        calphas_modes.buildHessian(calphas)
        calphas_modes.calcModes(n_modes=max_modes)
        modes = prody.extendModel(calphas_modes, calphas, molecule, norm=True)[0]
    else:
        modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        modes.buildHessian(molecule)
        modes.calcModes(n_modes=max_modes)
    return modes


def gaussian_modes(path):
    """
    Read the modes
    Create a prody.modes instance

    Parameters
    ----------
    path : str
        gaussian frequencies output path

    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    gaussian_parser = Gaussian(path).parse()
    shape = gaussian_parser.vibdisps.shape
    modes_vectors = gaussian_parser.vibdisps.reshape(shape[0], shape[1]*shape[2]).T
    modes_frequencies = np.abs(gaussian_parser.vibfreqs)
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


def chunker(end, n):
    """
    divide end integers in closed groups of n
    """
    for i in range(0, end-n+1, n):
        yield i+1, i+n
    if end % n:
        yield end-end % n+1, end


GROUPERS = {
    'residues': group_by_residues,
    'mass': group_by_mass,
    'calpha': 'calpha',
    '': None
}
