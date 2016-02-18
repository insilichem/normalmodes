from nmodMMTKinter import nmodMMTKinter
import numpy
class conf_change(nmodMMTKinter):
	def __init__(self, mols, molfile,option):
		self.molecules = mols
		self.molfile = molfile
		self.option = option
		self._mol1 = nmodMMTKinter(self.molecules, callback = self._loadMMTKCoordinates)

	def _loadMMTKCoordinates(self, analysis):
		self._mol1.loadMMTKCoordinates()
		self._mol2 = nmodMMTKinter(self.molfile, callback = self._loadMMTKCoordinates2)

	def _loadMMTKCoordinates2(self, analysis2):
		self._mol2.loadMMTKCoordinates()
		self.analysis()

	def analysis(self):

		from MMTK import Units
		from MMTK.NormalModes import EnergeticModes, VibrationalModes

		self.conf_mol1 = self._mol1.universe.copyConfiguration()
		tr, rms = self._mol1.universe.findTransformation(self.conf_mol1)
		self._mol1.universe.applyTransformation(tr)
		self.conf_mol2 = self._mol1.universe.copyConfiguration()
		self._mol1.universe.setConfiguration(self.conf_mol1)
		if self.option == "Using Vibrational Modes":
			modes = VibrationalModes(self._mol1.universe, 300.*Units.K)
		else:
			modes = EnergeticModes(self._mol1.universe, 300.*Units.K)

		diff = (self.conf_mol2 - self.conf_mol1).scaledToNorm(1.)

		plot = ModePlot(modes, diff)
		plot2 = ModePlot2(modes,diff, self.option)

from plotdialog import PlotDialog
class ModePlot(PlotDialog):

	def __init__(self, modes, diff):
		PlotDialog.__init__(self)
		self.mode_numbers = numpy.arange(6,len(modes))

		self.overlaps = [modes.rawMode(i).dotProduct(diff)**2
				for i in self.mode_numbers]

		self.subplot = self.add_subplot(1,1,1)
		self._displayData()

	def _displayData(self):

		ax = self.subplot
		ax.clear()
		lines = ax.plot(self.mode_numbers, numpy.add.accumulate(self.overlaps),
				c='r', marker='o',mfc='b')
		ax.set_xlabel("Mode Number")
		ax.set_ylabel("Cumulative squared overlap")
		ax.set_title("Cumulative sum overlapping of the modes")
		ax.grid(True)
		self.draw()

class ModePlot2(PlotDialog):

	def __init__(self, modes, diff, option):
		PlotDialog.__init__(self)
		self.modes = modes
		self.option = option
		self.mode_numbers = numpy.arange(6,len(self.modes))

		self.overlaps = [self.modes.rawMode(i).dotProduct(diff)**2
				for i in self.mode_numbers]

		self.subplot = self.add_subplot(1,1,1)
		if self.option == "Using Vibrational Modes":
			self.registerPickHandler(self._onPick)
		else:
			self.registerPickHandler(self._onPick2)
		self._displayData()

	def _displayData(self):

		ax = self.subplot
		ax.clear()
		lines = ax.plot(self.mode_numbers[:100], self.overlaps[:100],
				c='r', picker = True, marker='o', mfc='b')
		ax.set_xlabel("Mode number")
		ax.set_ylabel("Squared overlap")
		ax.set_title("Contribution of the first 100 nmodes")
		ax.grid(True)
		self.draw()

	def _onPick(self, event):
		from OnPicker import Picker
		Picker(self.modes[int(event.ind[0])+6])
		print event.ind[0]

	def _onPick2(self, event):
		print event.ind[0]

