# This class is needed to interface with MovieDialog
class NormalModeTraj:
	def __len__(self):
		return len(self.molecule.coordSets)
	def __getitem__(self, key):
		return None

from chimera.baseDialog import ModelessDialog
class Picker(ModelessDialog):
	"""Display MMTK normal modes by reconstructing a Chimera
	(trajectory) model from the MMTK universe and one particular
	normal mode."""

	buttons = ("Show", "Quit")
	title = "Normal Modes"
	oneshot = True

	def __init__(self, modes, *args, **kw):
		# We can override window title with
		# self.title = XXX
		# if we only knew what XXX should be.
		# Create proxy normal mode objects to simplify use
		# Save modes.  Currently we only use it to get
		# at the universe.
		self.mode = modes
		self.sf = float("1.0")
		self.movieDialog = None
		ModelessDialog.__init__(self, *args, **kw)
		self._createMovieDialog()

	def fillInUI(self, parent):
		import Pmw
		self.scaleWidget = Pmw.Counter(parent,
					labelpos = 'w',
					label_text = 'Scale factor:',
					label_justify = 'right',
					entryfield_value = '1.0',
					datatype = {'counter':'real'},
					entryfield_validate = {'validator':'real',
					'min' : '1.0', 'max' : '100.0'},
					increment=1.0)

		self.scaleWidget.pack(expand=False, fill="x")

	def Show(self):

		self.sf = float(self.scaleWidget.get())
		self._updateTrajectory()

	def Quit(self):

		if self.movieDialog:
			self.movieDialog.Quit()
			self.movieDialog = None
		self.Close()


	def _createMovieDialog(self):
		# First create the reference molecule
		from MMTK2Molecule import convert
		m, self._mmtk2chimera = convert(self.mode.universe,
								defaultCS=1)

		# Next we need to add the extra coordinate sets showing
		# normal mode extreme positions.  The trajectory consists
		# of 4 coordinate sets.	 CS 1 was created above and is
		# the reference coordinate set.	 CS 2 and 4 are the two
		# extreme positions while CS 3 is another reference.
		# When played back, you can see the periodic motion.
		# Initially, we set CS 2 and CS 4 to be reference also.
		cs2 = m.newCoordSet(2)			# + displacement
		cs3 = m.newCoordSet(3)			# reference
		cs4 = m.newCoordSet(4)			# - displacement
		for a in m.atoms:
			p = a.coord()
			a.setCoord(p, cs2)
			a.setCoord(p, cs3)
			a.setCoord(p, cs4)

		# Finally we need to create the MovieDialog for the
		# trajectory we just made
		ensemble = NormalModeTraj()
		ensemble.name = "Normal mode"
		keys = m.coordSets.keys()
		ensemble.startFrame = min(keys)
		ensemble.endFrame = max(keys)
		ensemble.molecule = m
		from Movie.gui import MovieDialog
		self.movieDialog = MovieDialog(ensemble)
		self._updateTrajectory()

	def _updateTrajectory(self):
		# Assume that the trajectory has been created and has
		# exactly 4 coordinate sets.  See above.  We update
		# only the two non-reference coordinate sets.
		from chimera import Vector
		mode = self.mode
		m = self.movieDialog.ensemble.molecule
		cs1 = m.findCoordSet(1)			# reference
		cs2 = m.findCoordSet(2)			# + displacement
		cs4 = m.findCoordSet(4)			# - displacement
		for ma in self.mode.universe.atomList():
			x, y, z = mode[ma] * 10 * float(self.sf)
			v = Vector(x, y, z)
			a = self._mmtk2chimera[ma]
			p = a.coord(cs1)
			a.setCoord(p + v, cs2)
			a.setCoord(p - v, cs4)

