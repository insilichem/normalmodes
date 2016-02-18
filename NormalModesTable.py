import chimera
import prody
from nmodMMTKinter import nmodMMTKinter
#import prody
#import Savenmod

from chimera.baseDialog import ModelessDialog
from chimera import Vector
from chimera import UserError
from chimera import replyobj
from chimera.extension import manager
from SimpleSession import SAVE_SESSION
from SimpleSession import sessionID
from SimpleSession.save import pickled
#from SimpleSession import idLookup
from SimpleSession import END_RESTORE_SESSION
from CGLtk.Table import SortableTable
import Pmw
import Tkinter
import base
from Movie.gui import MovieDialog
#from Savenmod.gui import SaveNMDialog
from StringIO import StringIO

# This class is only needed for the GUI.  The SortableTable
# API is much simpler if column values are attributes (like
# "frequency", "index" and "active")


class _NormalModeProxy(object):

    def __init__(self, index, active, freq, displacements, dialog):
        self._dialog = dialog
        self._active = active
#       self.active = active
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

#   buttons = ( "Show", "Quit", "Save", "Move" )
#   Show happens automatically when user selects a mode
#   Save is handled by session-saving
#   Move is incompletely implemented
    buttons = ("Quit", )
    title = "Normal Modes"
    oneshot = True

    def __init__(self, mol=None, modes=None, sesData=None, *args, **kw):
        self.showVectors = False
        self.vectorModel = None
        self.sessionHandler = None
        self.movieDialog = None
        self.modeNumber = False
        self.__recursion = False
        self.counter = 0
        self.cyl = 0
        self.con = 0
        self.rel = 0
        self.vecSf = 0
        self.atoms = list()
        self.modeData = list()
        self.title = None
        self.molecule = None

        if modes:
            if not isinstance(modes, prody.RTB) and not isinstance(modes, prody.ANM):
                # Create proxy normal mode objects to simplify use
                # with SortableTable.
                print modes
                self.molecule = mol
                self.title = "Normal Modes for %s" % mol.name
                print self.title
                atomMap = dict()
                for ma in modes.universe.atomList():
                    a = ma.getAtomProperty(ma, "chimera_atom")
                    atomMap[a] = ma
                    self.atoms.append(a)
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
                print modes
                self.molecule = mol
                self.title = "Normal Modes for %s" % mol.name
                cs1 = mol.newCoordSet(0)
                for atom in mol.atoms:
                    atom.setCoord(atom.coord(), cs1)
                    self.atoms.append(atom)

                for index, mode in enumerate(modes):
                    active = False
                    displacements = list()
                    v = mode.getEigvec()
                    freq = str(mode.getEigval())
                    for i in xrange(0, len(v), 3):
                        x, y, z = v[i:i+3] * 10
                        displacements.append(Vector(x, y, z))
                    self.modeData.append(_NormalModeProxy(index,
                                                          active, freq,
                                                          displacements,
                                                          self))
                self._createMovieDialog(self.molecule)
        else:
            self._sessionRestore(sesData)
        ModelessDialog.__init__(self, *args, **kw)
        self.sessionHandler = chimera.triggers.addHandler(SAVE_SESSION,
                                                          self._sessionSaveCB, None)

    def fillInUI(self, parent):
        # scaleWidget allows user to scale the displacement
        # vector for display purposes.

        # addColumn arguments are: (column title, attribute associated
        # with column,  layout format).  format=None means use default
        # display string.  format=bool means display using checkbutton
        # that can be used to change the value of the variable.
        t = SortableTable(parent)
        t.pack(expand=True, fill="both")
        t.addColumn("Mode", "index", format="%d")
        t.addColumn("Frequency", "frequency", format=None)
        t.addColumn("Active", "active", format=bool)
        t.setData(self.modeData)
        t.launch()
        self.modesTable = t
        #self.hydro = False
        #self.resid = False
        self.representation = 1

#       self.scaleWidget = Tkinter.Scale(parent,
#                               label="Scale Factor",
#                               from_=1.0, to=50.0, orient="horizontal", resolution = 1.0,
#                               command = self.__scaleCB)
#       self.scaleWidget.set(0.1)
#       self.scaleWidget.pack(expand=True, fill="both")

        self.scaleWidget = Pmw.Counter(parent,
                                       entryfield_command=self.__scaleCB,
                                       autorepeat=True,
                                       entryfield_value=1,
                                       entryfield_validate=self._validate)

        self.scaleWidget.pack(expand=True, fill="both")

        self.optionVector = Pmw.RadioSelect(parent,
                                            buttontype='checkbutton',
                                            labelpos='w',
                                            label_text='Show displacement Vectors:',
                                            command=self._show)
        self.optionVector.add("")
        self.optionVector.pack(expand=False, fill='both')
        self.vectorGroup = Pmw.Group(parent, tag_text='Displacement Vectors')
        self.vectorGroup.pack(expand=False, fill='both')
        self.vectorGroup.pack_forget()

        # test for drawing vectors
        self.__vecTypeDef = Tkinter.StringVar()
        self.__vecTypeDef.set('all atoms')
        self.vecType = Pmw.OptionMenu(self.vectorGroup.interior(),
                                      labelpos='w',
                                      label_text='Representation :',
                                      menubutton_textvariable=self.__vecTypeDef,
                                      items=('all atoms', 'no hydrogens', 'by residues'),
                                      command=self.__vecType)
        self.vecType.pack(anchor='w', fill="both")

        self.vectorType = Pmw.RadioSelect(self.vectorGroup.interior(),
                                          buttontype='radiobutton',
                                          orient='vertical',
                                          command=self.__scaleCB2)
        self.vectorType.add("Positive Displacement")
        self.vectorType.add("Negative Displacement")
        self.vectorType.setvalue("Positive Displacement")
        self.vectorType.pack(expand=False, fill="both")

    # defining the arrows
        self.sf = Pmw.ScrolledFrame(self.vectorGroup.interior(),
                                    labelpos='n', label_text='Modifying Arrows',
                                    usehullsize=1,
                                    hull_width=400,
                                    hull_height=220)

        self.arrow = Pmw.RadioSelect(self.vectorGroup.interior(),
                                     buttontype='checkbutton',
                                     labelpos='w',
                                     label_text='Configuration arrows',
                                     command=self._arrowShow)
        self.arrow.add("")
        self.arrow.pack(expand=False, fill='both')

        # function per entrar infos sobre les amplituds de les flexes
        self.__vectorCyl = Pmw.Counter(self.sf.interior(),
                                       entryfield_command=self.__cylCB,
                                       autorepeat=True,
                                       label_text="Radius of the cylinder",
                                       labelpos="w",
                                       entryfield_value=0.03,
                                       datatype={'counter': 'real', 'separator': '.'},
                                       increment=0.01,
                                       entryfield_validate=self._validateCyl)
        self.__vectorCyl.pack(expand=False, fill="both")

        self.__vectorCon = Pmw.Counter(self.sf.interior(),
                                       entryfield_command=self.__conCB,
                                       autorepeat=True,
                                       label_text="Radius of the base of the cone",
                                       labelpos="w",
                                       entryfield_value=0.3,
                                       datatype={'counter': 'real', 'separator': '.'},
                                       increment=0.01,
                                       entryfield_validate=self._validateCon)
        self.__vectorCon.pack(expand=False, fill="both")

        self.__vectorRel = Pmw.Counter(self.sf.interior(),
                                       entryfield_command=self.__relCB,
                                       autorepeat=True,
                                       label_text="Relative size of the cylinder",
                                       labelpos="w",
                                       entryfield_value=0.8,
                                       datatype={'counter': 'real', 'separator': '.'},
                                       increment=0.01,
                                       entryfield_validate=self._validateRel)
        self.__vectorRel.pack(expand=False, fill="both")

        self.__var = Tkinter.StringVar()
        self.__var.set('medium purple')
        self.__vectorColor = Pmw.OptionMenu(self.sf.interior(),
                                            labelpos='w',
                                            label_text='Choose Color:',
                                            menubutton_textvariable=self.__var,
                                            items=(
                                                'red', 'blue', 'cyan', 'medium purple', 'yellow', 'grey', 'green'),
                                            command=self.__colorCB,
                                            )
        self.__vectorColor.pack(anchor='w', fill="both")

        self.__vectorScale = Pmw.Counter(self.sf.interior(),
                                         entryfield_command=self.__scaleCB,
                                         autorepeat=True,
                                         label_text="Vector Scale Factor",
                                         labelpos="w",
                                         entryfield_value=50,
                                         datatype={'counter': 'numeric'},
                                         increment=1,
                                         entryfield_validate=self._validateSf)
        self.__vectorScale.pack(expand=False, fill="both")

        self.saveOptionVector = Pmw.RadioSelect(parent,
                                                buttontype='checkbutton',
                                                labelpos='w',
                                                label_text='Save Modes in Session:',
                                                command=self._showsave)
        self.saveOptionVector.add("")
        self.saveOptionVector.pack(expand=False, fill='both')
        self.vecSaveGroup = Pmw.Group(parent, tag_text='Modes to save in session')
        self.vecSaveGroup.pack(expand=False, fill='both')

        self.modeText = Pmw.EntryField(self.vecSaveGroup.interior(),
                                       label_text=" modes",
                                       labelpos='e',
                                       validate={"validator": "real"})
        self.modeText.pack(expand=False, fill="both")
        self.vecSaveGroup.pack_forget()

    def __scaleCB(self):
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

    def __scaleCB2(self, val):
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

    def _validate(self, text):
        try:
            int(text)
        except:
            return -1
        if int(text) <= 0:
            return -1
        if self.counter == 1:
            self.counter -= 1
            self.__scaleCB()
        self.counter += 1
        return 1

    def _validateCyl(self, text):
        try:
            float(text)
        except:
            return -1

        if float(text) <= 0.00:
            return -1

        if self.cyl == 1:
            self.cyl -= 1
            self.__cylCB()
        self.cyl += 1
        return 1

    def _validateCon(self, text):
        try:
            float(text)
        except:
            return -1

        if float(text) <= 0.00:
            return -1

        if self.con == 1:
            self.con -= 1
            self.__conCB()
        self.con += 1
        return 1

    def _validateRel(self, text):
        try:
            float(text)
        except:
            return -1

        if float(text) <= 0.00:
            return -1
        if float(text) >= 0.99:
            return -1

        if self.rel == 1:
            self.rel -= 1
            self.__relCB()
        self.rel += 1
        return 1

    def _validateSf(self, text):
        try:
            int(text)
        except:
            return -1

        if int(text) <= 0:
            return -1
        if float(text) >= 600:
            return -1

        if self.vecSf == 1:
            self.vecSf -= 1
            self.__scaleCB()
        self.vecSf += 1
        return 1

    def __vecType(self, tag):
        which = None
        for nmp in self.modeData:
            if not nmp.active:
                continue
            if which is not None:
                raise UserError("Only one mode may be selected")
            else:
                which = nmp

        if self.vecType.getcurselection() == 'all atoms':
            print "Plot all the atoms arrows"
            self.representation = 1
        if self.vecType.getcurselection() == 'no hydrogens':
            print "Plot heavy atoms arrows"
            self.representation = 2
        if self.vecType.getcurselection() == 'by residues':
            print "Plot mass weighted arrow for each residue"
            print "origin of the arrow is the center of mass of the residue"
            self.representation = 3

        if self.showVectors:
            self._drawVector(which)
        self._updateTrajectory(which, None, None)

    def __cylCB(self):
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

    def __conCB(self):
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

    def __relCB(self):
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

    def _modeNumber(self, tag, state):
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
            self.sf.pack(padx=5, pady=3, fill='both', expand=1)
        else:
            self.sf.pack_forget()

    def _showsave(self, tag, state):
        if state:
            self.vecSaveGroup.pack(expand=False, fill="both")
        else:
            self.vecSaveGroup.pack_forget()

    def _show(self, tag, state):
        if state:
            self.showVectors = True
            self.Show()
            # self.vectorGroup.expand()
            self.vectorGroup.pack(expand=False, fill="both")
        else:
            self.showVectors = False
            chimera.openModels.close(self.__bildModel)
            # self.vectorGroup.collapse()
            self.vectorGroup.pack_forget()

    def Show(self):
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
            # self._drawVector(which,hydro)
        self._updateTrajectory(which, None, None)

    def Quit(self):
        if self.sessionHandler:
            chimera.triggers.deleteHandler(SAVE_SESSION,
                                           self.sessionHandler)
            self.sessionHandler = None
        if self.movieDialog:
            self.movieDialog.Quit()
            self.movieDialog = None
        self.destroy()

    def Move(self):
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
        cs2 = m.findCoordSet(1)         # + displacement
        if cs2 is None:
            cs2 = m.newCoordSet(1)
        cs3 = m.findCoordSet(2)         # + displacement
        if cs3 is None:
            cs3 = m.newCoordSet(2)
        cs4 = m.findCoordSet(3)         # + displacement
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
        self.movieDialog = MovieDialog(ensemble, externalEnsemble=True)

    def _updateTrajectory(self, nmp, triggerName, ignore):
        # Assume that the trajectory has been created and has
        # exactly 4 coordinate sets.  See above.  We update
        # only the two non-reference coordinate sets.

        if self.__recursion:
            return
        self.__recursion = True

        m = self.movieDialog.ensemble.molecule
        cs1 = m.findCoordSet(0)         # reference
        cs2 = m.findCoordSet(1)         # + displacement
        cs4 = m.findCoordSet(3)         # - displacement

        sf = int(self.scaleWidget.get())

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

        bild = StringIO()  # Define a temporary file to write the vector
        # Put the color of the bild arrows as an initial info
        bild.write(".color %s\n" % self.__vectorColor.getcurselection())
        sf = int(self.__vectorScale.get())  # get the scaling factor of the arrow

        if self.representation == 1:
            for ma in m.atoms:
                x = dict_modes[ma][0] * sf
                y = dict_modes[ma][1] * sf
                z = dict_modes[ma][2] * sf
                v = Vector(x, y, z)
                if str(v) == "0 0 0":
                    pass
                else:
                    p = ma.coord(cs1)
                    if self.vectorType.getcurselection() == "Positive Displacement":
                        print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(
                            p[2]+v[2])+" %.2f %.2f %.2f" % (float(self.__vectorCyl.get()), float(self.__vectorCon.get()), float(self.__vectorRel.get()))
                    else:
                        print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]-v[0])+" "+str(p[1]-v[1])+" "+str(
                            p[2]-v[2])+" %.2f %.2f %.2f" % (float(self.__vectorCyl.get()), float(self.__vectorCon.get()), float(self.__vectorRel.get()))

        if self.representation == 2:
            for ma in m.atoms:
                if str(ma.element) == "H":
                    v = Vector(0, 0, 0)
                else:
                    x = dict_modes[ma][0] * sf
                    y = dict_modes[ma][1] * sf
                    z = dict_modes[ma][2] * sf
                    v = Vector(x, y, z)
                    if str(v) == "0 0 0":
                        pass
                    else:
                        p = ma.coord(cs1)
                        if self.vectorType.getcurselection() == "Positive Displacement":
                            print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(
                                p[2]+v[2])+" %.2f %.2f %.2f" % (float(self.__vectorCyl.get()), float(self.__vectorCon.get()), float(self.__vectorRel.get()))
                        else:
                            print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]-v[0])+" "+str(p[1]-v[1])+" "+str(
                                p[2]-v[2])+" %.2f %.2f %.2f" % (float(self.__vectorCyl.get()), float(self.__vectorCon.get()), float(self.__vectorRel.get()))

        if self.representation == 3:
            for residues in m.residues:
                v = Vector(0., 0., 0.)
                p = residues.labelCoord()  # Es el punt initial del vector.
                # Com a prova he posat que sigui la posicio del label.
                # Esta a prop del cm pero s'ha de canviar per el centre de massa exacte.
                for atoms in residues.atoms:
                    x = dict_modes[atoms][0] * sf/10
                    y = dict_modes[atoms][1] * sf/10
                    z = dict_modes[atoms][2] * sf/10
                    v = v+Vector(x, y, z)
                if self.vectorType.getcurselection() == "Positive Displacement":
                    print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]+v[0])+" "+str(p[1]+v[1])+" "+str(
                        p[2]+v[2])+" %.2f %.2f %.2f" % (float(self.__vectorCyl.get()), float(self.__vectorCon.get()), float(self.__vectorRel.get()))
                else:
                    print >> bild, ".arrow "+str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(p[0]-v[0])+" "+str(p[1]-v[1])+" "+str(
                        p[2]-v[2])+" %.2f %.2f %.2f" % (float(self.__vectorCyl.get()), float(self.__vectorCon.get()), float(self.__vectorRel.get()))

        bild.seek(0)
        self.__bildModel = chimera.openModels.open(bild, subid=10, type="Bild", identifyAs="Vectors",
                                                   temporary=True)
        # sameAs=m)
        bild.close()
        self.__recursion = False
    #
    # Code for saving and restoring session
    #

    def _sessionSaveCB(self, triggerName, myData, sessionFile):
        sesData = dict()
        sesData["version"] = 1
        sesData["molecule"] = sessionID(self.molecule)
        sesData["atoms"] = [sessionID(a) for a in self.atoms]
        modeData = list()
        if self.modeNumber == False:
            for nmp in self.modeData:
                modeData.append((nmp.active, nmp.frequency,
                                 [v.data() for v in nmp.displacements]))
        else:
            n = 0
            while n < int(self.modeText.getvalue()):
                print n, self.modeText.getvalue()
                modeData.append((self.modeData[n].active, self.modeData[n].frequency,
                                 [v.data() for v in self.modeData[n].displacements]))
                n += 1
        sesData["modeData"] = modeData
        print >> sessionFile, "nmtData = %s" % pickled(sesData)
        print >> sessionFile, """
try:
    from nmod import restoreNMTSession
    restoreNMTSession(nmtData)
except:
    reportRestoreError("Error restoring NormalModesTable interface")
"""

    def _sessionRestore(self, data):
        if not data:
            return
        version = data.get("version")
        try:
            f = self._restoreVersionMap[version]
        except KeyError:
            replyobj.error(
                "Version %d of NormalModesTable not supported "
                "in this version of Chimera.\n" % version)
        else:
            f(self, data)

    def _restoreNMTSession_v1(self, data):
        from SimpleSession import idLookup
        self.molecule = idLookup(data["molecule"])
        self.title = "Normal Modes for %s" % self.molecule.name
        self.atoms = [idLookup(sid) for sid in data["atoms"]]
        self.modeData = list()
        for index, (active, freq, dl) in enumerate(data["modeData"]):
            self.modeData.append(_NormalModeProxy(index,
                                                  active, freq,
                                                  [Vector(*v) for v in dl],
                                                  self))
        self._restoreHandler = chimera.triggers.addHandler(
            END_RESTORE_SESSION,
            self._restoreMovieDialog, None)

    def _restoreMovieDialog(self, triggerName, myData, triggerData):
        chimera.triggers.deleteHandler(END_RESTORE_SESSION,
                                       self._restoreHandler)
        del self._restoreHandler

        for inst in manager.instances:
            if not isinstance(inst, MovieDialog):
                continue
            if inst.ensemble.molecule is self.molecule:
                self.movieDialog = inst
                break

    _restoreVersionMap = {
        1:  _restoreNMTSession_v1,
    }


def restoreNMTSession(data):
    NormalModesTableDialog(sesData=data)
