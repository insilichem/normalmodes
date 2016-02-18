from chimera.baseDialog import ModelessDialog

# This class is only needed for the GUI.  The SortableTable
# API is much simpler if column values are attributes (like
# "frequency", "index" and "active")

class _NormalModeProxy(object):

	def __init__(self, index, active, freq, displacements, dialog):
		self._dialog = dialog
		self._active = active
#		self.active = active
		self.index = index
		self.frequency = freq
		self.displacements = displacements

	#
	# Hack to make only one mode selected at any time
	#
	def get_active(self):
		return self._active
	def set_active(self, a):
		self._active = a
		self._dialog._activeCB(self)
	active = property(get_active, set_active)

# This class is needed to interface with MovieDialog
class NormalModeTraj:
	def __len__(self):
		return len(self.molecule.coordSets)
	def __getitem__(self, key):
		return None

class NormalModesTableDialog(ModelessDialog):
	"""Display MMTK normal modes by reconstructing a Chimera
	(trajectory) model from the MMTK universe and one particular
	normal mode."""

#	buttons = ( "Show", "Quit", "Save", "Move" )
#	Show happens automatically when user selects a mode
#	Save is handled by session-saving
#	Move is incompletely implemented
	buttons = ( "Quit", )
	title = "Normal Modes"
	oneshot = True

	def __init__(self, mol=None, modes=None, sesData=None, *args, **kw):
		self.showVectors = False
		self.vectorModel = None
		self.sessionHandler = None
		self.movieDialog = None
		self.modeNumber = False
		self.__recursion = False
		if modes:
			# Create proxy normal mode objects to simplify use
			# with SortableTable.
			self.molecule = mol
			self.title = "Normal Modes for %s" % mol.name
			self.atoms = list()
			atomMap = dict()
			for ma in modes.universe.atomList():
				a = ma.getAtomProperty(ma, "chimera_atom")
				atomMap[a] = ma
				self.atoms.append(a)
			from chimera import Vector
			self.modeData = list()
			for index, mode in enumerate(modes):
				active = False
				freq = str(mode.frequency)
				displacements = list()
				for a in self.atoms:
					ma = atomMap[a]
					x, y, z = mode[ma] * 10
					displacements.append(Vector(x, y, z))
				self.modeData.append(_NormalModeProxy(index,
								active, freq,
								displacements,
								self))
			self._createMovieDialog(self.molecule)
		else:
			self._sessionRestore(sesData)
		ModelessDialog.__init__(self, *args, **kw)
		from SimpleSession import SAVE_SESSION
		import chimera
		self.sessionHandler = chimera.triggers.addHandler(SAVE_SESSION,
						self._sessionSaveCB, None)

	def fillInUI(self, parent):
		from CGLtk.Table import SortableTable
		import Pmw
		import Tkinter
		# scaleWidget allows user to scale the displacement
		# vector for display purposes.
		self.scaleWidget = Tkinter.Scale(parent,
		                      	label="Vector Scale Factor",
		                      	from_=1.0, to=50.0, orient="horizontal", resolution = 1.0,
		                      	command = self.__scaleCB)
		self.scaleWidget.set(0.1)
		self.scaleWidget.pack(expand=True, fill="both")
		t = SortableTable(parent)
		t.pack(expand=True, fill="both")
		# addColumn arguments are: (column title, attribute associated
		# with column,  layout format).  format=None means use default
		# display string.  format=bool means display using checkbutton
		# that can be used to change the value of the variable.
		t.addColumn("Mode", "index", format="%d")
		t.addColumn("Frequency", "frequency", format=None)
		t.addColumn("Active", "active", format=bool)
		t.setData(self.modeData)
		t.launch()
		self.modesTable = t
		self.hydro = False
		self.resid = False

		self.optionVector = Pmw.RadioSelect(parent,
						buttontype = 'checkbutton',
						labelpos='w',
						label_text = 'Show displacement Vectors:',
						command = self._show)
		self.optionVector.add("")
		self.optionVector.pack(expand=False, fill='both')
		self.vectorGroup = Pmw.Group(parent, tag_text='Displacement Vectors')
		self.vectorGroup.pack(expand=False, fill='both')
       		self.vectorGroup.pack_forget()
		# Part of the interface to define 
		# if hydrogen atoms will have there vectors depicted
		self.hydroType = Pmw.RadioSelect(self.vectorGroup.interior(),
						buttontype = 'checkbutton',
						labelpos='w',
						label_text = 'Hydrogen displacements:',
						command = self.__hydroType)
		self.hydroType.add("")
		self.hydroType.pack(expand = False, fill = "both")
                # test for drawing vectors
		self.__vectorColor= Pmw.OptionMenu(self.vectorGroup.interior(), 
						labelpos='w', 
						label_text='type of representation :', 
						menubutton_textvariable=self.__var,
						items = ('all atoms','all atoms - no hydrogens','center of mass'),
						command =self.__vecType,
						)
		self.__vecType.pack(anchor='w', fill="both")

		# Part of the interface to define 
		# if vectors are described as weighted on residues
		self.residType = Pmw.RadioSelect(self.vectorGroup.interior(),
						buttontype = 'checkbutton',
						labelpos='w',
						label_text = 'Residue displacements:',
						command = self.__residType)
		self.residType.add("")
		self.residType.pack(expand = False, fill = "both")

		self.vectorType = Pmw.RadioSelect(self.vectorGroup.interior(),
						buttontype = 'radiobutton',
						orient = 'vertical',
						command = self.__scaleCB)
		self.vectorType.add("Positive Displacement")
		self.vectorType.add("Negative Displacement")
		self.vectorType.setvalue("Positive Displacement")
		self.vectorType.pack(expand = False, fill = "both")

                # defining the arrows
		self.sf = Pmw.ScrolledFrame(self.vectorGroup.interior(),
                labelpos = 'n', label_text = 'Modifying Arrows',
                usehullsize = 1,
                hull_width = 400,
                hull_height = 220,
        	)

		self.arrow = Pmw.RadioSelect(self.vectorGroup.interior(),
						buttontype = 'checkbutton',
						labelpos = 'w',
						label_text = 'modifying arrows',
						command = self._arrowShow) 
		self.arrow.add("")
		self.arrow.pack(expand=False, fill='both')
		####


		#function per entrar infos sobre les amplituds de les flexes
		self.__vectorCyl=Tkinter.Scale(self.sf.interior(),
						label="Radius of the cylinder", 
						from_=0.03, 
						to=1, 
						orient="horizontal", 
						resolution=0.01, 
						command = self.__cylCB)
		self.__vectorCyl.set(0.1)
		self.__vectorCyl.pack(expand=False, fill="both")

		self.__vectorCon=Tkinter.Scale(self.sf.interior(),label="Radius of the base of the cone", from_=0.03, to=1, orient="horizontal", resolution=0.01, command = self.__conCB)
		self.__vectorCon.set(0.3)
		self.__vectorCon.pack(expand=False, fill="both")

		self.__vectorRel=Tkinter.Scale(self.sf.interior(),label="Relative size of the cylinder ", from_=0.1, to=0.9, orient="horizontal", resolution=0.1, command = self.__relCB)
		self.__vectorRel.set(0.8)
		self.__vectorRel.pack(expand=False, fill="both")

		self.__var = Tkinter.StringVar()
     		self.__var.set('medium purple')
		self.__vectorColor= Pmw.OptionMenu(self.sf.interior(), 
						labelpos='w', 
						label_text='Choose Color:', 
						menubutton_textvariable=self.__var,
						items = ('red','blue','cyan','medium purple','yellow','grey','green'),
						command =self.__colorCB,
						)
		self.__vectorColor.pack(anchor='w', fill="both")
	
		self.__vectorScale = Tkinter.Scale(self.sf.interior(),
						label="Vector Scale Factor",
						from_=0.1, to=500.0, orient="horizontal", resolution = 0.1,
						command = self.__scaleCB)
		self.__vectorScale.set(200)
		self.__vectorScale.pack(expand=False, fill="both")
		
		self.saveOptionVector = Pmw.RadioSelect(parent,
						buttontype = 'checkbutton',
						labelpos='w',
						label_text = 'Save Modes in Session:',
						command = self._showsave)
		self.saveOptionVector.add("")
		self.saveOptionVector.pack(expand=False, fill='both')
		self.vecSaveGroup = Pmw.Group(parent, tag_text='Modes to save in session')
		self.vecSaveGroup.pack(expand=False, fill='both')
#		self.entryModeNumber = Pmw.RadioSelect(self.vecSaveGroup.interior(),
#						label_text="Save only the first ",
#						labelpos='w',
###						buttontype = 'checkbutton',
						#command = self._modeNumber)
		#self.entryModeNumber.add("")
		#self.entryModeNumber.pack(expand=False,fill="x")
		self.modeText = Pmw.EntryField(self.vecSaveGroup.interior(), 
						label_text=" modes",
						labelpos = 'e',
						validate={"validator":"real"})
		self.modeText.pack(expand=False,fill="both")
		self.vecSaveGroup.pack_forget()

	def __scaleCB(self, val):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)

# 03/06/2013 - hydroType: function that defines whether to display or not hydrogen vectors

	def __hydroType(self,tag,state):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp
		if state:
			self.hydro = True
			print "Hydrogen vectors will be depicted"
		else:
			self.hydro = False
			print "Hydrogen vectors will be hidden"

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)

	def __residType(self,tag,state):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp
		if state:
			self.resid = True
			print "The vectors will be weighted on the residues"
			if self.hydro == True:
				self.hydroType.invoke("")
		else:
			self.resid = False
			print "The vectors will not be weighted on the residues"

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)


	def __cylCB(self, val):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)

	def __conCB(self, val):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)
	def __relCB(self, val):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)

	def __colorCB(self, color):
		from chimera import UserError
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp

		if self.showVectors:
			self._drawVector(which)
		self._updateTrajectory(which, None, None)
	
	def _modeNumber(self,tag,state):
		if state:
			self.modeNumber = True
		else:
			self.modeNumber = False
	
	def _activeCB(self, selected):
		#
		# Hack to make sure only one mode is selected at any time
		#
		if selected.active:
			changed = False
			for nmp in self.modeData:
				if nmp is not selected and nmp.active:
					nmp.active = False
					changed = True
			if changed:
				self.modesTable.refresh()
		self.Show()

	def _arrowShow(self, tag, state):
		if state:
			self.sf.pack(padx = 5, pady = 3, fill = 'both', expand = 1)
		else:
			self.sf.pack_forget()

	def _showsave(self, tag, state):
		if state:
			self.vecSaveGroup.pack(expand=False,fill="both")
		else:
			self.vecSaveGroup.pack_forget()

	def _show(self, tag, state):
		if state:
			self.showVectors = True
			self.Show()
			#self.vectorGroup.expand()
			self.vectorGroup.pack(expand=False, fill="both")
		else:
			self.showVectors = False
			import chimera
			chimera.openModels.close(self.__bildModel)
			#self.vectorGroup.collapse()
			self.vectorGroup.pack_forget()
	def Show(self):
		from chimera import UserError
		# For now, only allow one mode to be selected
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp
		sf = float(self.scaleWidget.get())
		self._updateTrajectory(which, None, None)
		if self.showVectors:
			self._drawVector(which)
			#self._drawVector(which,hydro)
		self._updateTrajectory(which, None, None)

	def Quit(self):
		if self.sessionHandler:
			from SimpleSession import SAVE_SESSION
			import chimera
			chimera.triggers.deleteHandler(SAVE_SESSION,
							self.sessionHandler)
			self.sessionHandler = None
		if self.movieDialog:
			self.movieDialog.Quit()
			self.movieDialog = None
		self.destroy()

	def Move(self):
		from chimera import UserError
		# For now, only allow one mode to be selected
		which = None
		for nmp in self.modeData:
			if not nmp.active:
				continue
			if which is not None:
				raise UserError("Only one mode may be selected")
			else:
				which = nmp
		if which is None:
			raise UserError("Please select normal mode to show.")
		sf = float(self.scaleWidget.get())
		self._moveTrajectory(which, sf)

	def _createMovieDialog(self, m):
		# We need to add the extra coordinate sets showing
		# normal mode extreme positions.  The trajectory consists
		# of 4 coordinate sets.  CS 1 was created above and is
		# the reference coordinate set.  CS 2 and 4 are the two
		# extreme positions while CS 3 is another reference.
		# When played back, you can see the periodic motion.
		# Initially, we set CS 2 and CS 4 to be reference also.
		cs2 = m.findCoordSet(1)			# + displacement
		if cs2 is None:
			cs2 = m.newCoordSet(1)
		cs3 = m.findCoordSet(2)			# + displacement
		if cs3 is None:
			cs3 = m.newCoordSet(2)
		cs4 = m.findCoordSet(3)			# + displacement
		if cs4 is None:
			cs4 = m.newCoordSet(3)
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
		ensemble.startFrame = min(keys)+1
		ensemble.endFrame = max(keys)+1
		ensemble.molecule = m
		from Movie.gui import MovieDialog
		self.movieDialog = MovieDialog(ensemble, externalEnsemble=True)

	def _updateTrajectory(self, nmp, triggerName, ignore):
		# Assume that the trajectory has been created and has
		# exactly 4 coordinate sets.  See above.  We update
		# only the two non-reference coordinate sets.
		
		if self.__recursion:
			return
		self.__recursion = True


		m = self.movieDialog.ensemble.molecule
		cs1 = m.findCoordSet(0)			# reference
		cs2 = m.findCoordSet(1)			# + displacement
		cs4 = m.findCoordSet(3)			# - displacement

		sf = self.scaleWidget.get()
	
		if nmp is not None:
			for a, d in zip(self.atoms, nmp.displacements):
				v = d * sf
				p = a.coord(cs1)
				a.setCoord(p + v, cs2)
				a.setCoord(p - v, cs4)
		else:
			for a in self.atoms:
				p = a.coord(cs1)
				a.setCoord(p, cs2)
				a.setCoord(p, cs4)

		self.__recursion = False

	def Save(self):
		import Savenmod
		from Savenmod.gui import SaveNMDialog
		SaveNMDialog(self.modes)

	#
	# Extension manager stuff
	#
	def emHide(self):
		self.Cancel()

	def emName(self):
		return self.title

	def emQuit(self):
		self.Quit()

	def emRaise(self):
		self.enter()

	def _drawVector(self, nmp):
		from chimera import Vector
		import chimera
		
		if self.__recursion:
			return
		else:
			self.__recursion = True
		try:
			chimera.openModels.close(self.__bildModel)
		except:
			pass
		
		m = self.movieDialog.ensemble.molecule
		cs1 = m.findCoordSet(1)
		mode = nmp.displacements
		dict_modes = dict()
		for ma, ma1 in zip(mode, m.atoms):
			dict_modes[ma1] = ma
		
		from StringIO import StringIO
	
		bild = StringIO() # Define a temporary file to write the vector
		bild.write(".color %s\n" % self.__vectorColor.getcurselection()) # Put the color of the bild arrows as an initial info
		sf = self.__vectorScale.get() # get the scaling factor of the arrow

		if self.resid==True:
			bild2 = StringIO() # Define a temporary file to write the vector
			bild2.write(".color %s\n" % self.__vectorColor.getcurselection()) # Put the color of the bild arrows as an initial info
			for residues in m.residues:
				v=Vector(0.,0.,0.)
				for atoms in residues.atoms:
					x = dict_modes[atoms][0] * sf
					y = dict_modes[atoms][1] * sf
					z = dict_modes[atoms][2] * sf
					v=v+Vector(x,y,z)
				p=residues.labelCoord() # Es el punt initial del vector. Com a prova he posat que sigui la posicio del label. Esta a prop del cm pero s'ha de canviar per el centre de massa exacte. 
				print >> bild2, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(p[2]+v[2])+" %.2f %.2f %.2f" % (self.__vectorCyl.get(),self.__vectorCon.get(),self.__vectorRel.get())
			bild2.seek(0)
			self.__bildModel = chimera.openModels.open(bild2, type="Bild", identifyAs="Vectors residue",
							temporary = True, sameAs=m)
			bild2.close()
			self.__recursion = False
		else:
			

		if self.hydro==False:
			for ma in m.atoms:
				if str(ma.element)=="H":
					v = Vector(0,0,0)
				else:
					x = dict_modes[ma][0] * sf
					y = dict_modes[ma][1] * sf
					z = dict_modes[ma][2] * sf
					v = Vector(x,y,z)
				#TEST JD PER IMPRIMIR NOMES VECTOR NO HYDROGEN			
				# hauria de ser cap a 0 no 0 just"
					if str(v) == "0 0 0":
						pass
					else:
						p = ma.coord(cs1)
						if self.vectorType.getcurselection() == "Positive Displacement":
							print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(p[2]+v[2])+" %.2f %.2f %.2f" % (self.__vectorCyl.get(),self.__vectorCon.get(),self.__vectorRel.get())
						else:
							print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]-v[0])+" "+str(p[1]-v[1])+" "+str(p[2]-v[2])+" %.2f %.2f %.2f" % (self.__vectorCyl.get(),self.__vectorCon.get(),self.__vectorRel.get())
			bild.seek(0)
			self.__bildModel = chimera.openModels.open(bild, type="Bild", identifyAs="Vectors",
								temporary = True, sameAs=m)
			bild.close()
			self.__recursion = False

		else:
			for ma in m.atoms:
				x = dict_modes[ma][0] * sf
				y = dict_modes[ma][1] * sf
				z = dict_modes[ma][2] * sf
				v = Vector(x,y,z)
				#TEST JD PER IMPRIMIR NOMES VECTOR NO HYDROGEN			
				# hauria de ser cap a 0 no 0 just"
				if str(v) == "0 0 0":
					pass
				else:
					p = ma.coord(cs1)
					if self.vectorType.getcurselection() == "Positive Displacement":
						#print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(p[2]+v[2])+" 0.03 0.08 0.8"
						print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(p[2]+v[2])+" %.2f %.2f %.2f" % (self.__vectorCyl.get(),self.__vectorCon.get(),self.__vectorRel.get())
					else:
						print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]-v[0])+" "+str(p[1]-v[1])+" "+str(p[2]-v[2])+" %.2f %.2f %.2f" % (self.__vectorCyl.get(),self.__vectorCon.get(),self.__vectorRel.get())
			bild.seek(0)
			self.__bildModel = chimera.openModels.open(bild, type="Bild", identifyAs="Vectors",
								temporary = True, sameAs=m)
			bild.close()
			self.__recursion = False

	#
	# Code for saving and restoring session
	#
	def _sessionSaveCB(self, triggerName, myData, sessionFile):
		from SimpleSession import sessionID
		sesData = dict()
		sesData["version"] = 1
		sesData["molecule"] = sessionID(self.molecule)
		sesData["atoms"] = [ sessionID(a) for a in self.atoms ]
		modeData = list()
		if self.modeNumber == False:
			for nmp in self.modeData:
				modeData.append((nmp.active, nmp.frequency,
							[ v.data() for v in nmp.displacements ]))
		else:
			n=0
			while n < int(self.modeText.getvalue()):
				print n,self.modeText.getvalue()
				modeData.append((self.modeData[n].active, self.modeData[n].frequency,
							[v.data() for v in self.modeData[n].displacements]))
				n+=1
		sesData["modeData"] = modeData
		from SimpleSession.save import pickled
		print >> sessionFile, "nmtData = %s" % pickled(sesData)
		print >> sessionFile, """
try:
	from nmod import restoreNMTSession
	restoreNMTSession(nmtData)
except:
	reportRestoreError("Error restoring NormalModesTable interface")
"""

	def _sessionRestore(self, data):
		version = data["version"]
		try:
			f = self._restoreVersionMap[version]
		except KeyError:
			from chimera import replyobj
			replyobj.error(
				"Version %d of NormalModesTable not supported "
				"in this version of Chimera.\n" % version)
		else:
			f(self, data)

	def _restoreNMTSession_v1(self, data):
		from SimpleSession import idLookup
		self.molecule = idLookup(data["molecule"])
		self.title = "Normal Modes for %s" % self.molecule.name
		self.atoms = [ idLookup(sid) for sid in data["atoms"] ]
		self.modeData = list()
		from chimera import Vector
		for index, (active, freq, dl) in enumerate(data["modeData"]):
			self.modeData.append(_NormalModeProxy(index,
						active, freq,
						[ Vector(*v) for v in dl ],
						self))
		import chimera
		from SimpleSession import END_RESTORE_SESSION
		self._restoreHandler = chimera.triggers.addHandler(
						END_RESTORE_SESSION,
						self._restoreMovieDialog, None)

	def _restoreMovieDialog(self, triggerName, myData, triggerData):
		import chimera
		from SimpleSession import END_RESTORE_SESSION
		chimera.triggers.deleteHandler(END_RESTORE_SESSION,
						self._restoreHandler)
		del self._restoreHandler

		from chimera.extension import manager
		from Movie.gui import MovieDialog
		for inst in manager.instances:
			if not isinstance(inst, MovieDialog):
				continue
			if inst.ensemble.molecule is self.molecule:
				self.movieDialog = inst
				break

	_restoreVersionMap = {
		1:	_restoreNMTSession_v1,
	}

def restoreNMTSession(data):
	NormalModesTableDialog(sesData=data)
