from nmodMMTKinter import nmodMMTKinter
import CalcNormalModes as cnm
import CoarseGrainAlgorithm as CGAlg


class nmod(nmodMMTKinter):

    def __init__(self, mols, proc,
                 prodyalgorithm=None, n_algorithm=None):
        self.molecules = mols

        if prodyalgorithm == (None or 'Residues'):
            algorithm = CGAlg.alg1
        elif prodyalgorithm == 'Mas':
            algorithm = CGAlg.alg2
        self.modes = cnm.calc_normal_modes(self.molecules[0], algorithm, n_algorithm)
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
        NormalModesTableDialog(self.molecules[0], self.modes)

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
