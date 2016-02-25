# --- UCSF Chimera Copyright ---
# Copyright (c) 2006 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---

import chimera
from chimera.baseDialog import ModelessDialog
import Pmw
import Tkinter
from DockPrep.gui import DockPrepDialog
from Tkinter import *
from chimera.widgets import MoleculeScrolledListBox
import chimera.dialogs
from MMMD.gui import mmmdDialog


class NMDialog(ModelessDialog):
    name = "Normal Modes Calculation"
    buttons = ("Run", "Close")

    def fillInUI(self, parent):
        self.lenJon = False
        self.cutOff = False
        self.var_subSpace = False
        self.proc = "std"
        parent.columnconfigure(1, weight=1)
        row = 0

        self.BackGround = None
        self.fixAtoms = None
        self.subspace = None

        #
        # Calculation procedure menu
        #

        self.NMProcedureMenu = Pmw.OptionMenu(parent,
                                              labelpos='w',
                                              label_text="Force Field:",
                                              items=[
                                                  'Amber', 'Elastic Network', 'Ca - not yet available', 'Prody'],
                                              command=self._FFmenu)
        self.NMProcedureMenu.grid(column=0, row=row, sticky='w')

        self.AlgorithmMenu = Pmw.OptionMenu(parent,
                                            labelpos='w',
                                            label_text='Algorithm:',
                                            items=['Residues', 'Mas'])
        # command=self._algorithm_menu)

        self.AlgorithmDialog = Pmw.EntryField(parent,
                                              validate={'validator': 'numeric'},
                                              labelpos='w',
                                              label_text='n:',
                                              entry_width=5)
        row += 1

        #
        # Molecule Minimization
        #
        self.MPLabel = Tkinter.Label(parent,
                                     text="Minimize structure:")
        self.MPLabel.grid(column=0, row=row, sticky='w')
        self.MinimizerButton = Tkinter.Button(parent,
                                              text="Proceed",
                                              command=lambda: chimera.dialogs.display(
                                                  mmmdDialog.name))
        self.MinimizerButton.grid(column=1, row=row, sticky='w')
        row += 1

        #
        # Model selection
        #
        self.molList = MoleculeScrolledListBox(parent,
                                               labelpos='w',
                                               label_text="Select model:   ",
                                               listbox_selectmode="extended")
        self.molList.grid(column=0, row=row, columnspan=3,
                          sticky='nsew')
        parent.rowconfigure(row, weight=1)
        row += 1

        #
        # DockPrep
        #
        self.dockPrepGroup = Pmw.Group(parent,
                                       tag_text="System Setup")
        self.dockPrepGroup.grid(column=0, row=row, columnspan=3,
                                sticky='nsew')
        parent.rowconfigure(row, weight=1)
        row += 1

        self.dockButton = Tkinter.Button(self.dockPrepGroup.interior(), text="Run",
                                         command=lambda: chimera.dialogs.display(DockPrepDialog.name))
        self.memoryType = Pmw.RadioSelect(self.dockPrepGroup.interior(), orient="vertical",
                                          buttontype='radiobutton', pady=0)
        self.memoryType.add("set", text="Memorize options chosen in"
                            " subsequent dialogs", anchor="w", justify="left")
        self.memoryType.add("use", text="Use previously memorized options,"
                            " if any", anchor="w", justify="left")
        self.memoryType.add("none", text="Neither memorize nor use memorized"
                            " options", anchor="w", justify="left")
        self.memoryType.grid(row=1, column=0)
        from DockPrep.prefs import prefs, MEMORIZED_SETTINGS
        if "Minimize" in prefs[MEMORIZED_SETTINGS]:
            self.memoryType.setvalue("use")
        else:
            self.memoryType.setvalue("set")
        self.dockButton.grid(column=0, row=0, sticky='nsew')

        #
        # Procedure Options
        #
        self.ProcOptions = Pmw.Group(parent, tag_text="Flexibility")
        self.ProcOptions.grid(column=0, row=row, columnspan=3, sticky='nsew')

        self.var1 = Tkinter.StringVar()
        self.var1.set('All Atoms')
        self.procMenu = Pmw.OptionMenu(self.ProcOptions.interior(),
                                       labelpos='w',
                                       label_text="Choose: ",
                                       menubutton_textvariable=self.var1,
                                       items=[
                                           'All Atoms', 'Rigid Body Residues', 'Fixed Selection'],
                                       command=self._procOptions)
        self.procMenu.grid(row=0, column=0, sticky='nsew')
        self.var2 = Tkinter.StringVar()
        self.var2.set('-')

        self.fixMenu = Pmw.OptionMenu(self.ProcOptions.interior(),
                                      labelpos='w',
                                      label_text="Flexible Space:",
                                      menubutton_textvariable=self.var2,
                                      items=['-'])
        self.fixMenu.setitems(["selected", "unselected"], index="selected")
        self.fixMenu.grid(row=1, column=0, sticky='nsew')
        self.fixSetButton = Tkinter.Button(self.ProcOptions.interior(),
                                           text="Set", command=self._setAtoms)
        self.fixSetButton.grid(column=1, row=1, sticky='w')
        self.fixSetButton.config(state=DISABLED)
        row += 1

        #
        # Force Field Options
        #
        self.ffOptions = Pmw.Group(parent, tag_text="Forcefield Options")
        self.ffOptions.grid(column=0, row=row, columnspan=3, sticky='nsew')
        self.cutOffOption = Pmw.RadioSelect(self.ffOptions.interior(),
                                            buttontype='checkbutton',
                                            command=self._cutOff)
        self.cutOffOption.add('')
        self.cutOffOption.grid(column=0, row=0, sticky='nsew')
        self.cutOffDialog = Pmw.EntryField(self.ffOptions.interior(),
                                           validate={'validator': 'real'},
                                           labelpos='w',
                                           label_text='Electrostatic Options:')
        self.cutOffDialog.grid(row=0, column=1, sticky='nsew')
        self.var1 = Tkinter.StringVar()
        self.var1.set('direct')
        self.cutOffMenu = Pmw.OptionMenu(self.ffOptions.interior(),
                                         menubutton_textvariable=self.var1,
                                         items=['direct', 'cutoff', 'ewald', 'screened', 'multipole'])
        self.cutOffMenu.grid(row=0, column=2, sticky='nsew')
        self.var2 = Tkinter.StringVar()
        self.var2.set('direct')

        self.lenOption = Pmw.RadioSelect(self.ffOptions.interior(),
                                         buttontype='checkbutton',
                                         command=self._lenJon)
        self.lenOption.add('')
        self.lenOption.grid(column=0, row=1, sticky='nsew')
        self.lenDialog = Pmw.EntryField(self.ffOptions.interior(),
                                        validate={'validator': 'real'},
                                        labelpos='w',
                                        label_text='Lennar-Jones Options:')
        self.lenDialog.grid(row=1, column=1, sticky='nsew')
        self.lenMenu = Pmw.OptionMenu(self.ffOptions.interior(),
                                      menubutton_textvariable=self.var2,
                                      items=['direct', 'cutoff'])
        self.lenMenu.grid(row=1, column=2, sticky='nsew')
        row += 1

        #
        # Conformational mode analysis
        #
        self.groupAnalysis = Pmw.Group(parent,
                                       tag_text="Conformational Mode Analysis")
        self.groupAnalysis.grid(column=0, row=row, columnspan=3,
                                sticky='nsew')
        parent.rowconfigure(row, weight=1)
        row += 1
        self.analysisCheck = Pmw.RadioSelect(
            self.groupAnalysis.interior(),
            buttontype='checkbutton',
            command=self._analysisCheckCB)
        self.analysisCheck.add("Perform Conformational Analysis")
        self.analysisCheck.grid(column=0, row=0, sticky='nsew')
        self.modeAnalysis = Pmw.RadioSelect(
            self.groupAnalysis.interior(),
            buttontype='radiobutton')
        self.modeAnalysis.add("Using Energetic Modes")
        self.modeAnalysis.add("Using Vibrational Modes")
        self.modeAnalysis.setvalue("Using Energetic Modes")
        self.modeAnalysis.grid(column=0, row=1, sticky="nsew")
        self.molList2 = MoleculeScrolledListBox(
            self.groupAnalysis.interior(),
            labelpos='w',
            label_text="Select model\n"
            " to compare:",
            listbox_selectmode="extended")
        self.molList2.grid(column=0, row=2, sticky='nsew')
        self.analysisButton = Tkinter.Button(
            self.groupAnalysis.interior(),
            text="Analysis",
            command=self.analysis,
            state=DISABLED)
        self.analysisButton.grid(column=0, row=3, sticky='nsew')

        #
        # Team information
        #
        self.teamName = Pmw.Group(parent)
        self.teamName.grid(column=0, row=row, columnspan=3,
                           sticky="nsew")
        self.teamNameInfo = Tkinter.Label(self.teamName.interior(),
                                          text="Interface designed by V. Munoz-Robles "
                                          "and J.-D.Marechal "
                                          "and... Jordi Guasp\n"
                                          "InsiliChem\n")
        self.teamNameInfo.grid(column=0, row=0, sticky="ew")
        row += 1

    def _FFmenu(self, tag):
        if tag == "Elastic Network":
            self.ffOptions.grid_remove()
            self.AlgorithmMenu.grid_remove()
            self.AlgorithmDialog.grid_remove()
        elif tag == 'Prody':
            self.ffOptions.grid_remove()
            self.AlgorithmMenu.grid(row=0, column=1, sticky='w')
            self.AlgorithmDialog.grid(row=0, column=2, sticky='ew', padx=5)
        else:
            self.ffOptions.grid()
            self.AlgorithmMenu.grid_remove()
            self.AlgorithmDialog.grid_remove()

    # def _algorithm_menu(self, tag):
    #     if tag == 'Residues':
    #         pass
    #     elif tag == 'Mas':
    #         pass
    #     else:
    #         pass

    def _cutOff(self, tag, state):
        if state:
            self.cutOff = True
        else:
            self.cutOff = False

    def _lenJon(self, tag, state):
        if state:
            self.lenJon = True
        else:
            self.lenJon = False

    def _setAtoms(self):
        from chimera import runCommand as run
        if self.fixMenu.getvalue() == 'selected':
            run("sel invert sel")
            self.fixAtoms = chimera.selection.currentAtoms()
        elif self.fixMenu.getvalue() == 'unselected':
            self.fixAtoms = chimera.selection.currentAtoms()
        elif self.fixMenu.getvalue() == 'Backbone':
            run("sel @n,ca,c,o,h,ha")
            # run("sel invert sel")
            self.fixAtoms = chimera.selection.currentAtoms()
        elif self.fixMenu.getvalue() == 'Backbone-minimal':
            run("sel @n,ca,c,h,ha")
            # run("sel invert sel")
            self.fixAtoms = chimera.selection.currentAtoms()
        elif self.fixMenu.getvalue() == 'Calpha':
            run("sel @ca")
            # run("sel invert sel")
        self.fixAtoms = chimera.selection.currentAtoms()
        run("color blue ~sel")
        run("namesel fixAtoms")
        run("~sel")

    def _procOptions(self, opt):
        if opt == "Fixed Selection":
            self.fixSetButton.config(state=NORMAL)
        elif opt == "Rigid Body Residues":
            self.fixSetButton.config(state=DISABLED)
            self.fixAtoms = "subspace"
        else:
            self.fixSetButton.config(state=DISABLED)
            self.fixAtoms = None

#   def _runningOptionsCB(self,tag):
#
#       if self.RunningOptions.getcurselection()=="Foreground":
#           self.OutputNameGenerate.config(state=DISABLED)
#           self.BackGround = None
#           from chimera import help
#           help.register(self.RunningOptions, balloon="Foreground option will held chimera until the caculation is complete")
#       else:
#           self.OutputNameGenerate.config(state=NORMAL)
#           self.BackGround = True

    def _analysisCheckCB(self, tag, state):

        if state:
            self.analysisButton.config(state=NORMAL)
            # self.modeAnalysis.config(state = NORMAL)
        else:
            self.analysisButton.config(state=DISABLED)
            # self.modeAnalysis.config(state = DISABLED)

    def analysis(self):
        option = self.modeAnalysis.getcurselection()
        from confChange import conf_change
        analysis = conf_change(self.molList.getvalue(),
                               self.molList2.getvalue(), option)

    def generate(self):

        self.filename = self.OutputNameLabel.getvalue()
        from chimera import UserError
        if not self.filename:
            raise UserError("Please enter name for output file")

        raise UserError("Not implemented yet")

    def Run(self):
        esOptions = None
        ljOptions = None
        if self.cutOff:
            value = float(self.cutOffDialog.get())/10
            esOptions = {self.cutOffMenu.getcurselection(): value}
        if self.lenJon:
            value = float(self.lenDialog.get())/10
            ljOptions = {self.lenMenu.getcurselection(): value}
        self.molecules = self.molList.getvalue()

        # if self.fixAtoms == None and self.proc != "std":
        # from chimera import UserError
        # raise UserError("No atoms were selected")

        if not self.molecules:
            from chimera import UserError
            raise UserError("No molecules selected")

        from base import nmod
#       if self.BackGround:
#           filename = self.OutputNameLabel.get()
#       else:
#           filename = None
        filename = None
        prodyalgorithm = None
        n_algorithm = None

        if self.NMProcedureMenu.getcurselection() == "Amber":
            runOpt = "ffm"
        elif self.NMProcedureMenu.getcurselection() == "Elastic Network":
            runOpt = "enm"
        elif self.NMProcedureMenu.getcurselection() == "Prody":
            runOpt = 'prd'
            prodyalgorithm = self.AlgorithmMenu.getcurselection()
            n_algorithm = int(self.AlgorithmDialog.get())

        self.mi = nmod(self.molecules, runOpt, self.proc, filename,
                       self.BackGround, fix=self.fixAtoms, memorize=self.memoryType.getvalue(),
                       ljOptions=ljOptions, esOptions=esOptions,
                       prodyalgorithm=prodyalgorithm, n_algorithm=n_algorithm)


from chimera import dialogs
dialogs.register(NMDialog.name, NMDialog)
