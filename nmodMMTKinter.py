from MMMD.MMTKinter import *

class nmodMMTKinter(MMTKinter):

	def __init__(self, mols, *args, **kw):

		MMTKinter.__init__(self, mols, *args, **kw)

	def _makeUniverse(self):
		MMTKinter._makeUniverse(self)

	def NM_ffm(self, proc, filename, BackGround = None, fix=None):
		# Novel interface working on the 25/09/2013
		import MMTK
		from MMTK.NormalModes import VibrationalModes
		from MMTK.Trajectory import Trajectory
		from MMTK.Subspace import Subspace
		if not fix:
			modes = VibrationalModes(self.universe)
			return modes
		elif fix == "subspace":
			import chimera
			from MMTK.Subspace import RigidMotionSubspace
			from MMTK.NormalModes import SubspaceNormalModes
			from MMTK.Collections import Collection

			def createCollection(atoms,subspace):
				col=Collection()
				for atom in atoms:
					col.addChemicalObject(atom)
				return col
			
			space=list()
			for res in self.mols[0].residues:
					space.append(createCollection([self.atomMap[atom] for atom in res.atoms],subspace))

			rigid=RigidMotionSubspace(self.universe,space)
			modes=SubspaceNormalModes(self.universe,rigid)
			return modes

		else:
			from MMTK.NormalModes import NormalModes, SubspaceNormalModes
			from MMTK.Subspace import RigidMotionSubspace
			from MMTK import Subspace
			import chimera
			selected_atoms = [self.atomMap[a] for a in fix]
			self.subspace = RigidMotionSubspace(self.universe,selected_atoms)
			modes = SubspaceNormalModes(self.universe, self.subspace)
			return modes

	def modesVisualization(self, filename, universe = None):
		import MMTK
		from nmod import NormalModesTable
		modes = MMTK.load("%s.nmod" % filename)
		if universe:
			modes.universe = universe
		else:
			modes.universe = self.universe
		NormalModesTable.NormalModesTableDialog(self.mols[0], modes)

	def generate_nmod(self, filename, runningOptions):

		import MMTK
		saveData = self.saveData()
		self.deleteChimeraAtoms()
		self.mark_atoms()

		if runningOptions == "FFM":
			self.writeScript_NM_ffm(filename)
			MMTK.save(self.universe, "%s.mmtk" % filename)
			self.restoreChimeraAtoms(saveData)
		if runningOptions == "MCM":
			import chimera
			if not chimera.selection.currentAtoms():
				from tkMessageBox import showerror
				showerror("Normal Modes MCM Calculation","No atoms were selected\nPlease select the atoms that meant to be fixed")
				pass
			else:
				selected_atoms = chimera.selection.invertCurrent()
				selected_atoms = [self.atomMap[a] for a in chimera.selection.currentAtoms()]
			for a in self.universe.atomList():
				if a in selected_atoms:
					a.mark = 1
				else:
					a.mark = 0
			MMTK.save(self.universe, "%s.mmtk" % filename)
			self.restoreChimeraAtoms(saveData)
			self.writeScript_NM_mcm(filename)
		if runningOptions == "ENM":
			from MMTK import InfiniteUniverse
			from MMTK.ForceFields import DeformationForceField

			universe = InfiniteUniverse(DeformationForceField())
			saveData = self.saveData()
			self.deleteChimeraAtoms()
			self.mark_atoms()
			for mm in self.molecules:
				self.universe.removeObject(mm)
				universe.addObject(mm)
			MMTK.save(self.universe, "%s.mmtk" % filename)
			self.restoreChimeraAtoms(saveData)
			self.writeScript_NM_enm(filename)
		if runningOptions == "C-alpha":
			MMTK.save(self.universe, "%s.mmtk" % filename)
			self.restoreChimeraAtoms(saveData)
			self.writeScript_NM_calpha(filename)

	def writeScript_NM_ffm(self, filename):
		input = open("%s.py" % filename, "w")
		print >> input, "import MMTK"
		print >> input, "from MMTK.NormalModes import NormalModes"
		print >> input, "universe = MMTK.load(\"%s.mmtk\")" % filename
		print >> input, "modes = NormalModes(universe)"
		print >> input, "MMTK.save(modes, \"%s.nmod\")" % filename
		print >> input, "MMTK.save(universe, \"%s.mmtk\")" % filename
		input.close()


	def NM_enm(self, proc, filename, BackGround = None, fix=None, subspace=None):

		from MMTK import InfiniteUniverse
		from MMTK.Subspace import Subspace
		from MMTK.Subspace import RigidMotionSubspace
		from MMTK.ForceFields import DeformationForceField
		from MMTK.FourierBasis import FourierBasis, estimateCutoff
		from MMTK.NormalModes import NormalModes, SubspaceNormalModes, VibrationalModes
		from MMTK.Collections import Collection

		# Construct system
		universe = InfiniteUniverse(DeformationForceField())

		self.setFixed("none")

		if fix == "subspace":

			def createCollection(atoms,subspace):
				col=Collection()
				for atom in atoms:
					col.addChemicalObject(atom)
				return col
			
			space=list()
			for res in self.mols[0].residues:
				space.append(createCollection([self.atomMap[atom] for atom in res.atoms],subspace))

		elif not fix:
			pass

		else:
			selected_atoms = [self.atomMap[a] for a in fix]

		for mm in self.molecules:
			self.universe.removeObject(mm)
			universe.addObject(mm)
			timestamp("Added model %s" % mm.name)

		# Find a reasonable basis set size and cutoff

		if not fix:
			nbasis = max(10, universe.numberOfAtoms()/5)
			cutoff, nbasis = estimateCutoff(universe, nbasis)
			print "Calculating %d low-frequency modes." % nbasis
			if cutoff is None:
			# Do full normal mode calculation
				modes = NormalModes(universe, sub)
			else:
			# Do subspace mode calculation with Fourier basis
				subspace = FourierBasis(universe, cutoff)
				modes = SubspaceNormalModes(universe, subspace)
			return modes
		if fix == "subspace":
			subspace=RigidMotionSubspace(universe,space)
			modes=SubspaceNormalModes(universe,subspace)
			return modes
		else:
			subspace = RigidMotionSubspace(universe,selected_atoms)
			modes = SubspaceNormalModes(universe, subspace)
			return modes			

	def writeScript_NM_enm(self, filename):
		input = open("%s.py" % filename, "w")
		print >> input, "import MMTK"
		print >> input, "from MMTK.NormalModes import NormalModes, SubspaceNormalModes"
		print >> input, "from MMTK.FourierBasis import FourierBasis, estimateCutoff"
		print >> input, "universe = MMTK.load(\"%s.mmtk\")" % filename
		print >> input, "nbasis = max(10, universe.numberOfAtoms()/5)"
		print >> input, "cutoff, nbasis = estimateCutoff(universe, nbasis)"
		print >> input, "if cutoff is None:"
		print >> input, "\tmodes = NormalModes(universe)"
		print >> input, "else:"
		print >> input, "\tsubspace = FourierBasis(universe, cutoff)"
		print >> input, "\tmodes = SubspaceNormalModes(universe, subspace)"
		print >> input, "MMTK.save(modes, \"%s.nmod\")" % filename
		print >> input, "MMTK.save(universe, \"%s.mmtk\")" % filename

		input.close()

	def NM_calpha(self):

		from MMTK import Units
		from MMTK import InfiniteUniverse
		from MMTK.ForceFields import CalphaForceField
		from MMTK.FourierBasis import FourierBasis, estimateCutoff
		from MMTK.NormalModes import EnergeticModes

		# Construct system

		universe = InfiniteUniverse(CalphaForceField(2.5))
		self.setFixed("none")

		for mm in self.molecules:
			self.universe.removeObject(mm)
			universe.addObject(mm)
			timestamp("Added model %s" % mm.name)

		# Find a reasonable basis set size and cutoff
		nbasis = max(10, universe.numberOfAtoms()/5)
		cutoff, nbasis = estimateCutoff(universe, nbasis)
		print "Calculating %d low-frequency modes." % nbasis


		if cutoff is None:
			# Do full normal mode calculation
			modes = EnergeticModes(universe,300.*Units.K)
		else:
			# Do subspace mode calculation with Fourier basis
			subspace = FourierBasis(universe, cutoff)
			modes = EnergeticModes(universe,300.*Units.K,subspace)

		return modes

	def deleteChimeraAtoms(self):

		for a in self.universe.atomList():
			del a.chimera_atom

		for a in self.molecules:
			del a.needParmchk

		for a in self.universe.objectList():
			del a.atomMap
			del a.chimeraMolecule

	def saveData(self):

		saveData = dict()
		for a in self.universe.objectList():
			saveData[a] = [a.atomMap, a.chimeraMolecule]
		return saveData

	def restoreChimeraAtoms(self, saveData):

		for ma in self.universe.atomList():
			for ca, ma in self.atomMap.iteritems():
				ma.chimera_atom = ca

		for a in self.universe.objectList():
			a.atomMap = saveData[a][0]
			a.chimeraMolecule = saveData[a][1]

	def changeUniverse(self):

		self.atomMap2 = dict()

		conf = self.universe.configuration()
		conf2 = self.universe2.configuration()

		for ma2, ma in zip(self.universe2.atomList(), self.universe.atomList()):
			ma2.chimera_atom = ma.chimera_atom
			self.atomMap2[ma2.chimera_atom] = ma2

		for m2, m in zip(self.universe2.objectList(), self.universe.objectList()):
			m2.chimeraMolecule = m.chimeraMolecule
			m2.atomMap = m.atomMap

		self.universe = self.universe2
		del self.universe2

		self.atomMap = self.atomMap2
		del self.atomMap2

		def update_md(conf):
			from chimera import Coord
			for  ma in self.universe.atomList():
				x, y, z = conf[ma] * 10
				ca = ma.getAtomProperty(ma, "chimera_atom")
				ca.setCoord(Coord(x, y, z))

		update_md(self.universe.configuration())

	def mark_atoms(self):
		n = 0
		for a in self.universe.atomList():
			a.id = n
			n += 1

	def newUniverse(self, filename):
		import MMTK
		self.universe2 = MMTK.load("%s.mmtk" % filename)
		self.changeUniverse()


