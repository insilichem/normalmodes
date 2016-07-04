import sys
import CalcNormalModes as cnm
import CoarseGrainAlgorithm as CGAlg
# SITE_PKGS = '/home/jguasp/.local/UCSF-Chimera64-1.10.2/lib/python2.7/site-packages'
# def remove_numpy_from_modules():
#     numpy_mods = [key for key in sys.modules.keys()
#                   if 'numpy' in key.lower()]
#     print 'Removing %d NumPy modules' % (len(numpy_mods),)
#     for numpy_mod in numpy_mods:
#         sys.modules.pop(numpy_mod)

# sys.path.insert(0, SITE_PKGS)
# remove_numpy_from_modules()
# import numpy
import Scientific
# import Numeric
# import numpy
# numpy.oldnumeric = Numeric
from nmodMMTKinter import nmodMMTKinter
# sys.path.pop(0)
# remove_numpy_from_modules()


class nmod(nmodMMTKinter):

    def __init__(self, mols, proc,
                 prodyalgorithm=None, n_algorithm=None,
                 LJ=False, mass_weighted=True):
        self.molecules = mols

        if prodyalgorithm == 'Residues':
            algorithm = CGAlg.alg1
        elif prodyalgorithm == 'Mas':
            algorithm = CGAlg.alg2
        else:
            algorithm = None

        self.modes = cnm.calc_normal_modes(self.molecules[0], algorithm, n_algorithm,
                                           LJ=LJ, mass_weighted=mass_weighted)
        # remove_numpy_from_modules()
        # sys.path.insert(0, SITE_PKGS)
        # import numpy

        self.visualization()

        self.modes = nmod

    # def needDockPrep(self):
    #     for m in self.molecules:
    #         if getattr(m, 'chargeModel', None) is None:
    #             return True
    #         for a in m.atoms:
    #             if getattr(a, 'charge', None) is None or getattr(a, 'gaffType', None) is None:
    #                 return True
    #     return False

    # def _loadCoord(self, nm):
    #     if isinstance(nm, nmodMMTKinter):
    #         nm.loadMMTKCoordinates()
    #     if self.runOpt == "ffm":
    #         self.modes = nm.NM_ffm(self.proc, self.filename, self.BackGround, fix=self.fix)
    #         if self.modes and not self.BackGround:
    #             self.visualization()
    #     elif self.runOpt == "enm":
    #         self.modes = nm.NM_enm(
    #             self.proc, self.filename, BackGround=self.BackGround, fix=self.fix)
    #         if not self.BackGround:
    #             self.visualization()
        # elif self.runOpt == "prd":
        #     # self.modes = calc modes with prody (self.proc, self.filename, BackGround
        #     # = self.BackGround, fix=self.fix)
        #     self.modes = cnm.calc_normal_modes(moldy)
        #     # modesprody to modes mmtk?
        #     if not self.BackGround:
        #         self.visualization()

    def visualization(self):

        from NormalModesTable import NormalModesTableDialog

        self.freq_dialog = NormalModesTableDialog(self.molecules[0], self.modes)

class gaussian_nmod(nmodMMTKinter):

    def __init__(self, mols, proc,
                 file):
        self.molecules = mols

        self.modes = gaussian2prody(file)

        self.visualization()

    def visualization(self):

        from NormalModesTable import NormalModesTableDialog
        NormalModesTableDialog(self.molecules[0], self.modes)

def gaussian2prody(file):
    pass

# _miCache = []


# def _find(molecules, exclres, nogui, addhyd, callback, memorize, cache, prep, esOptions, ljOptions):
#     global _miCache
#     # print "_find", id(_miCache), _miCache
#     mset = set(molecules)
#     for t in _miCache:
#         mols, exres, nm = t
#         if set(mols) == mset and exres == exclres:
#             callback(nm)
#             return

#     def cacheIt(nm, molecules=molecules, exclres=exclres,
#                 callback=callback, cache=_miCache):
#         # print "cacheIt", id(cache), cache
#         callback(nm)
#         cache.append((molecules, exclres, nm))
#         # print "end cacheIt", id(cache), cache
#     # print "creating MMTKinter instance", molecules
#     from nmodMMTKinter import nmodMMTKinter
#     nm = nmodMMTKinter(molecules, exclres=exclres,
#                        nogui=nogui, addhyd=addhyd,
#                        callback=(cacheIt if cache else callback),
#                        memorize=memorize, prep=prep, ljOptions=ljOptions, esOptions=esOptions)


# def _moleculeCheck(triggerName, data, mols):
#     # Remove all entries that refer to a molecule that is being closed
#     _removeFromCache(mols)


# def _removeFromCache(mols):
#     global _miCache
#     # print "Remove from cache", id(_miCache), _miCache
#     junk = []
#     for t in _miCache:
#         molecules = t[0]
#         for m in mols:
#             if m in molecules:
#                 junk.append(t)
#                 break
#     for t in junk:
#         _miCache.remove(t)
