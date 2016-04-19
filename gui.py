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

        self.AlgorithmMenu = Pmw.OptionMenu(parent,
                                            labelpos='w',
                                            label_text='Algorithm:',
                                            items=['Residues', 'Mas'])

        self.AlgorithmDialog = Pmw.EntryField(parent,
                                              validate={'validator': 'numeric'},
                                              labelpos='w',
                                              label_text='n:',
                                              entry_width=5)

        self.AlgorithmMenu.grid(row=row, column=0, sticky='w')
        self.AlgorithmDialog.grid(row=row, column=1, sticky='ew', padx=5)

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
        # Team information
        #
        self.teamName = Pmw.Group(parent)
        self.teamName.grid(column=0, row=row, columnspan=3,
                           sticky="nsew")
        self.teamNameInfo = Tkinter.Label(self.teamName.interior(),
                                          text="Interface designed by V. Munoz-Robles, "
                                          "J.-D.Marechal and... Jordi Guasp\n"
                                          "The computational Biotechnological "
                                          "Chemistry Team")
        self.teamNameInfo.grid(column=0, row=0, sticky="ew")
        row += 1

    def Run(self):
        self.molecules = self.molList.getvalue()

        if not self.molecules:
            from chimera import UserError
            raise UserError("No molecules selected")

        from base import nmod

        prodyalgorithm = None
        n_algorithm = None

        prodyalgorithm = self.AlgorithmMenu.getcurselection()
        n_algorithm = int(self.AlgorithmDialog.get())

        self.mi = nmod(self.molecules, self.proc,
                       prodyalgorithm=prodyalgorithm, n_algorithm=n_algorithm)


from chimera import dialogs
dialogs.register(NMDialog.name, NMDialog)
