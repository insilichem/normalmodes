#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import Tkinter as tk
from tkFileDialog import askopenfilename
import Pmw
import webbrowser as web
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MoleculeScrolledListBox
from libplume.ui import PlumeBaseDialog
from core import Controller
#from MMMD import mmmdDialog

ui = None
def showUI():
    if chimera.nogui:
        tk.Tk().withdraw()
    global ui
    if not ui:
        ui = NormalModesExtension()
    ui.enter()


class NormalModesExtension(PlumeBaseDialog):

    buttons = ('OK',)

    def __init__(self, *args, **kwargs):
        # GUI init
        self.title = 'Plume Normal Modes'
        self.modes_dialog = None

        # Variables
        self.var_input_choice = tk.StringVar()
        self.var_input_choice.set('prody')
        self.var_input_choice.trace('w', self._check_choice)

        # Fire up
        super(NormalModesExtension, self).__init__(with_logo=False, resizable=False,
                                                   *args, **kwargs)

    def fill_in_ui(self, parent):
        self.ui_input_frame = tk.LabelFrame(self.canvas, text='Select mode', padx=5, pady=5)
        self.ui_input_frame.pack(expand=True, fill='x', padx=5, pady=5)
        self.ui_input_choice_frame = tk.Frame(self.ui_input_frame)
        self.ui_input_choice_frame.grid(row=0)

        self.ui_input_choice_prody = tk.Radiobutton(self.ui_input_choice_frame,
                                                    variable=self.var_input_choice,
                                                    text='ProDy', value='prody')
        self.ui_input_choice_gaussian = tk.Radiobutton(self.ui_input_choice_frame,
                                                       variable=self.var_input_choice,
                                                       text='Gaussian', value='gaussian')
        self.ui_input_choice_prody.pack(side='left')
        self.ui_input_choice_gaussian.pack(side='left')
        self.ui_input_choice_prody.select()

    def Apply(self):
        if self.modes_dialog is None:
            self.modes_dialog = NormalModesConfigDialog(self,
                                    engine=self.var_input_choice.get())
        self.modes_dialog.enter()

    def OK(self):
        self.Apply()
        self.Close()

    def Close(self):
        global ui
        ui = None
        super(NormalModesExtension, self).Close()

    def _check_choice(self, *a, **kw):
        value = self.var_input_choice.get()
        if value == 'prody':
            try:
                import prody
            except ImportError as e:
                self.status('Prody is not properly installed!', color='red')
                print('Mode not available because {}'.format(e))
        elif value == 'gaussian':
            try:
                import cclib
            except ImportError as e:
                self.status('cclib must be installed!', color='red')
                print('Mode not available because {}'.format(e))


class NormalModesConfigDialog(PlumeBaseDialog):

    buttons = ('Run', 'Close')
    help = "https://github.com/insilichem/plume_normalmodes"
    VERSION = '0.0.1'
    VERSION_URL = "https://api.github.com/repos/insilichem/plume_normalmodes/releases/latest"


    def __init__(self, parent=None, engine='prody', *args, **kwargs):
        # GUI init
        self.title = 'Calculate Normal Modes with ' + engine.title()
        self.parent = parent
        self.engine = engine

        if engine == 'prody':
            self.fill_in_ui = self._fillInUI_Prody
        else:
            self.fill_in_ui = self._fillInUI_Gaussian

        self.lennard_jones = False
        self.mass_weighted = False

        super(NormalModesConfigDialog, self).__init__(resizable=False,
                                                      *args, **kwargs)

    def _fillInUI_Prody(self, parent):
        #
        # Algorithm selection
        #
        self.canvas.columnconfigure(1, weight=1)
        row = 0
        self.ui_algorithms_menu_lbl = tk.Label(self.canvas, text='Algorithm:',
                                               anchor='e')
        self.ui_algorithms_menu = Pmw.OptionMenu(self.canvas,
                                                 items=['Full atom',
                                                        'Extend from C-alpha',
                                                        'Group by residues',
                                                        'Group by mass',
                                                        'Group by graph'],
                                                 command=self._algorithms_menu_cb)

        self.ui_algorithms_param = Pmw.EntryField(self.canvas, entry_width=3,
                                                  validate={'validator': 'integer',
                                                            'min': 1},
                                                  labelpos='w', label_text='#:')
        self.ui_algorithms_param.setvalue(1)
        self.ui_algorithms_menu_lbl.grid(row=row, column=0, sticky='we', padx=3, pady=3)
        self.ui_algorithms_menu.grid(row=row, column=1, sticky='we', padx=3, pady=3, columnspan=2)


        #
        # Number of modes & Molecule Minimization
        #
        row += 1
        self.ui_n_modes_lbl = tk.Label(self.canvas, text='# Modes:', anchor='e')
        self.ui_n_modes = Pmw.EntryField(self.canvas, entry_width=3,
                                        validate={'validator': 'integer',
                                                  'min': 1})
        self.ui_n_modes.setvalue(20)
        self.ui_n_modes_lbl.grid(column=0, row=row, sticky='we', padx=3, pady=3)
        self.ui_n_modes.grid(column=1, row=row, sticky='news', padx=3)
        self.ui_minimize_btn = tk.Button(self.canvas, text="Minimize",
                                         )#command=lambda: chimera.dialogs.display(mmmdDialog.name))
        self.ui_minimize_btn.grid(column=2, row=row, sticky='we', padx=3, pady=3)

        #
        # Cutoff
        #
        row += 1
        self.ui_cutoff_lbl = tk.Label(self.canvas, text='Cutoff:', anchor='e')
        self.ui_cutoff = Pmw.EntryField(self.canvas, entry_width=3,
                                         validate={'validator': 'real',
                                                   'min': 4.0})
        self.ui_cutoff.setvalue(15.0)
        self.ui_cutoff_lbl.grid(column=0, row=row, sticky='we', padx=3, pady=3)
        self.ui_cutoff.grid(column=1, row=row, sticky='news', padx=3)

        #
        # Cutoff
        #
        row += 1
        self.ui_gamma_lbl = tk.Label(self.canvas, text='Gamma LJ:', anchor='e')
        self.ui_gamma = Pmw.EntryField(self.canvas, entry_width=3,
                                        validate={'validator': 'real',
                                                  'min': 0.0})
        self.ui_gamma.setvalue(1.0)
        self.ui_gamma_lbl.grid(column=0, row=row, sticky='we', padx=3, pady=3)
        self.ui_gamma.grid(column=1, row=row, sticky='news', padx=3)
        #
        # Optional Selections: Lennard-Jones and mass-weighted hessian
        #
        row += 1
        self.ui_extra_options = Pmw.Group(self.canvas, tag_text='Options')
        self.ui_extra_options.grid(column=0, row=row, columnspan=3, sticky='nsew',
                                   padx=5, pady=5)

        self.ui_extra_options_chk = Pmw.RadioSelect(self.ui_extra_options.interior(),
                                                    buttontype='checkbutton')
        self.ui_extra_options_chk.add('Mass-Weighted')
        self.ui_extra_options_chk.grid(column=0, row=0, sticky='we')

        #
        # Model selection
        #
        row += 1
        self.ui_molecules_lbl = tk.Label(self.canvas, text="Select model:",
                                         anchor='e')
        self.ui_molecules_lbl.grid(column=0, row=row, padx=3)
        self.ui_molecules = MoleculeScrolledListBox(self.canvas,
                                                    listbox_selectmode="single")
        self.ui_molecules.grid(column=1, row=row, columnspan=3, sticky='nsew',
                               padx=3, pady=3)

    def _fillInUI_Gaussian(self, parent):
        row = 0
        parent.columnconfigure(0, weight=1)
        self.ui_gaussian_grp = Pmw.Group(self.canvas,
                tag_text='Open Gaussian output file')
        self.ui_gaussian_grp.grid(column=0, row=row, columnspan=2,
                                  sticky='nsew', padx=5, pady=5)
        self.ui_gaussian_grp.columnconfigure(0, weight=1)
        self.ui_gaussian_grp.columnconfigure(1, weight=0)
        self.ui_gaussian_file_entry = Pmw.EntryField(self.ui_gaussian_grp.interior())
        self.ui_gaussian_file_entry.pack(side='left', expand=True,
                                         fill='both', padx=5, pady=5)

        self.ui_gaussian_btn = tk.Button(self.ui_gaussian_grp.interior(),
                                         text='...', command=self._load_file)
        self.ui_gaussian_btn.pack(side='right', padx=5, pady=5)

    def Apply(self):
        self.controller = Controller(self)
        self.vibrations = self.controller.run()

    def Run(self):
        self.Apply()
        self.Close()

    def _load_file(self, *a, **kw):
        path = askopenfilename()
        if path:
            self.ui_gaussian_file_entry.setvalue(path)

    def _algorithms_menu_cb(self, *a, **kw):
        value = self.ui_algorithms_menu.getvalue()
        if value.startswith('Group by'):
            self.ui_algorithms_param.grid(row=0, column=2, sticky='ew', padx=5)
            if value == 'Group by residues':
                self.ui_algorithms_param.configure(label_text='# residues:')
            elif value == 'Group by mass':
                self.ui_algorithms_param.configure(label_text='# groups:')
        else:
            self.ui_algorithms_param.grid_forget()


class NormalModesResultsDialog(PlumeBaseDialog):

    buttons = ('Close',)

    def __init__(self, parent=None, controller=None, *args, **kwargs):
        # GUI init
        self.title = 'Normal Modes Results'
        self.parent = parent
        self.controller = controller
        super(NormalModesResultsDialog, self).__init__(*args, **kwargs)

    def fill_in_ui(self, parent):
        pass

    def populate_data(self, frequencies=None):
        pass

    def plot_vectors(self, vectors=None):
        pass


class NormalModesMovieDialog(PlumeBaseDialog):

    buttons = ('Close')

    def __init__(self, parent=None, controller=None, *args, **kwargs):
        # GUI init
        self.title = 'Normal Modes Results'
        self.parent = parent
        self.controller = controller
        super(NormalModesMovieDialog, self).__init__(*args, **kwargs)

    def fill_in_ui(self, parent):
        pass
