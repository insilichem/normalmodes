#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import Tkinter as tk
from tkFileDialog import askopenfilename
import Pmw
import webbrowser as web
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MoleculeScrolledListBox
from plumesuite.ui import PlumeBaseDialog
from core import Controller
#from MMMD import mmmdDialog

ui = None
def showUI(callback=None):
    if chimera.nogui:
        tk.Tk().withdraw()
    global ui
    if not ui:
        ui = NormalModesExtension()
    ui.enter()
    if callback:
        ui.addCallback(callback)


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
        super(NormalModesExtension, self).__init__(self, *args, **kwargs)

    def fill_in_ui(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        self.ui_input_frame = tk.LabelFrame(self.canvas, text='Select mode', padx=5, pady=5)
        self.ui_input_frame.pack(expand=True, fill='x')
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
        """
        Default! Triggered action if you click on an Apply button
        """
        if self.modes_dialog is None:
            self.modes_dialog = NormalModesConfigDialog(self,
                                    engine=self.var_input_choice.get())
        self.modes_dialog.enter()

    def OK(self):
        """
        Default! Triggered action if you click on an OK button
        """
        self.Apply()
        self.Close()

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

        super(NormalModesConfigDialog, self).__init__(self, resizable=False,
                                                      *args, **kwargs)

    def _fillInUI_Prody(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        #
        # Algorithm selection
        #
        row = 0
        self.ui_algorithms_menu = Pmw.OptionMenu(self.canvas, labelpos='w',
                                                 label_text='Algorithm:',
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
        self.ui_algorithms_menu.grid(row=row, column=0, sticky='we', padx=3, pady=3)


        #
        # Number of modes & Molecule Minimization
        #
        row += 1
        self.ui_n_modes = Pmw.EntryField(self.canvas, entry_width=3,
                                        validate={'validator': 'integer',
                                                  'min': 1},
                                        labelpos='w', label_text='# Modes:')
        self.ui_n_modes.setvalue(20)
        self.ui_n_modes.grid(column=0, row=row, sticky='we', padx=3, pady=3)
        self.ui_minimize_btn = tk.Button(self.canvas, text="Minimize",
                                         )#command=lambda: chimera.dialogs.display(mmmdDialog.name))
        self.ui_minimize_btn.grid(column=1, row=row, sticky='we', padx=3, pady=3)


        # #
        # # Optional Selections: Lennard-Jones and mass-weighted hessian
        # #
        # row += 1
        # self.ui_extra_options = Pmw.Group(self.canvas, tag_text='Options')
        # self.ui_extra_options.grid(column=0, row=row, columnspan=2, sticky='nsew')

        # self.ui_extra_options_chk = Pmw.RadioSelect(self.ui_extra_options.interior(),
        #                                             buttontype='checkbutton')
        # self.ui_extra_options_chk.add('Lennard-Jones')
        # self.ui_extra_options_chk.add('Mass weighted hessian')
        # self.ui_extra_options_chk.grid(column=0, row=0, sticky='we')

        #
        # Model selection
        #
        row += 1
        self.ui_molecules = MoleculeScrolledListBox(self.canvas, labelpos='w',
                                                    label_text="Select model:",
                                                    listbox_selectmode="single")
        self.ui_molecules.grid(column=0, row=row, columnspan=2, sticky='nsew',
                               padx=3, pady=3)
        parent.rowconfigure(row, weight=1)

    def _fillInUI_Gaussian(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """

        row = 0
        self.ui_gaussian_grp = Pmw.Group(self.canvas,tag_text='Open Gaussian output file')
        self.ui_gaussian_grp.grid(column=0, row=row, columnspan=2,
                                  sticky='nsew', padx=5, pady=5)
        self.ui_gaussian_grp.columnconfigure(0, weight=1)
        self.ui_gaussian_grp.columnconfigure(1, weight=0)
        self.ui_gaussian_file_entry = Pmw.EntryField(self.ui_gaussian_grp.interior())
        self.ui_gaussian_file_entry.pack(side='left', expand=True, fill='x')

        self.ui_gaussian_btn = tk.Button(self.ui_gaussian_grp.interior(),
                                         text='...', command=self._load_file)
        self.ui_gaussian_btn.pack(side='right')

    def Apply(self):
        """
        Default! Triggered action if you click on an Apply button
        Change in core for apply_prody or apply_gaussian
        """
        self.controller = Controller(self)
        self.vibrations = self.controller.run()

    def Run(self):
        """
        Default! Triggered action if you click on an Run button
        """
        self.Apply()
        self.Close()

    def _load_file(self, *a, **kw):
        _hidden_files_fix()
        path = askopenfilename()
        if path:
            self.ui_gaussian_file_entry.setvalue(path)

    def _algorithms_menu_cb(self, *a, **kw):
        value = self.ui_algorithms_menu.getvalue()
        if value.startswith('Group by'):
            self.ui_algorithms_param.grid(row=0, column=1, sticky='ew', padx=5)
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
        super(NormalModesResultsDialog, self).__init__(self, *args, **kwargs)

    def fill_in_ui(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
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
        super(NormalModesMovieDialog, self).__init__(self, *args, **kwargs)

    def fill_in_ui(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        pass


class _HyperlinkManager:
    """
    from http://effbot.org/zone/tkinter-text-hyperlink.htm
    """

    def __init__(self, text):

        self.text = text

        self.text.tag_config("hyper", foreground="#367159", underline=1)

        self.text.tag_bind("hyper", "<Enter>", self._enter)
        self.text.tag_bind("hyper", "<Leave>", self._leave)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}

    def add(self, action):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        return "hyper", tag

    def _enter(self, event):
        self.text.config(cursor="hand2")

    def _leave(self, event):
        self.text.config(cursor="")

    def _click(self, event):
        for tag in self.text.tag_names(tk.CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag]()
                return
