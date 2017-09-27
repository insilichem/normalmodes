#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import Tkinter as tk
from tkFileDialog import askopenfilename
import Pmw
import webbrowser as web
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MoleculeScrolledListBox
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


class NormalModesExtension(ModelessDialog):
    
    default = None
    buttons = ('OK',)
    provideStatus = True
    statusPosition = 'left'
    help = 'https://www.insilichem.com'

    def __init__(self, *args, **kwarg):
        # GUI init
        self.title = 'Plume Normal Modes'
        self.modes_dialog = None

        # Variables
        self.var_input_choice = tk.StringVar()
        self.var_input_choice.set('prody')

        # Fire up
        ModelessDialog.__init__(self, resizable=False)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def fillInUI(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

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

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        chimera.extension.manager.deregisterInstance(self)
        self.destroy()
    
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


class NormalModesConfigDialog(ModelessDialog):

    buttons = ('Run', 'Close')
    default = None
    
    help = 'https://www.insilichem.com'

    def __init__(self, parent=None, engine='prody', *args, **kwarg):
        # GUI init
        self.title = 'Calculate Normal Modes with ' + engine.title()
        self.parent = parent
        self.engine = engine

        if engine == 'prody':
            self.fillInUI = self._fillInUI_Prody
        else:
            self.fillInUI = self._fillInUI_Gaussian

        self.lennard_jones = False
        self.mass_weighted = False

        # Fire up
        ModelessDialog.__init__(self, resizable=False)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def _team_logo(self, parent):
        parent = tk.Frame(parent)
        # InsiliChem copyright
        bg = chimera.tkgui.app.cget('bg')
        img_data = r"R0lGODlhZABrAOeeADhwWDBzWTlxWTF0WjpyWjtzWzx0XD11XT52Xj93X0B4YEF5YUJ6YkN7Yz98aEt6Y0R8ZEB9aUx7ZEV9ZU18ZUZ+Zk59Zkt+bE9+Z0x/bVB/aE2AblGAaU6Bb1KBak+CcFOCa1CDcVSDbFGEclWEbVKFc1OHdFSIdVuGdVWJdlyHdlaKd12Id1eLeF6JeFiMeV+KeWCLemGMe2KNfGOOfWSPfmWQf2aSgGeTgW6RgWiUgm+SgmmVg3CTg3GUhHKVhXOWhnSXh3WYiHaZiXebinici3mdjICcjXqejYGdjnufjoKej3ygj4OfkH2hkISgkYGhl4WhkoKimIaik4OjmYeklIWlm4illYmmloennYqnmIionouomYmpn4ypmoqqoI2qm4uroY6rnJGqopKro5OspJStpZWuppiuoJavp5ewqJixqZqyqpuzq5y0rJ21rZ62rqS1rp+3r6W2r6a3sKe4sai5sqm6tKq7tau8tqy9t62+uK6/ubW/uq/BurbAu7DCu7jBvLLDvLnCvbPEvbrDvrvEv7zFwL3Hwr7Iw7/JxMDKxcHLxsLMx8PNyMrMycTOycbPysvOys3Py87QzM/RztDSz9HT0NLU0dPV0tTW09XX1NbY1dja1v///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH5BAEKAP8ALAAAAABkAGsAQAj+ADsJHEiwoMGDCA1ienTITpgiNEJAMAAhU4YPK2wU8SKnECNKmxKKHEmyJME9BgCkaUQAgMuXMGPKBMBAU8uZOGfGoeIyTCaTQDtFUpBz5iMylIIqJRiJyqSURWMGAjonqssC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"
        img = tk.PhotoImage(data=img_data)
        logo = tk.Button(parent, image=img, background=bg, borderwidth=0,
                         activebackground=bg, highlightcolor=bg, cursor="hand2",
                         command=lambda *a: web.open_new(r"http://www.insilichem.com/"))
        logo.image = img
        text = tk.Text(parent, background=bg, borderwidth=0, height=4, width=30)
        hrefs = _HyperlinkManager(text)
        big = text.tag_config("big", font="-size 18", foreground="#367159")
        text.insert(tk.INSERT, "InsiliChem", "big")
        text.insert(tk.INSERT, "\nDeveloped by ")
        text.insert(tk.INSERT, "@jaimergp", hrefs.add(lambda *a: web.open_new(r"https://github.com/jaimergp")))
        text.insert(tk.INSERT, "\nat Maréchal Group, UAB, Spain")
        text.configure(state='disabled')
        logo.grid(row=5, column=0, sticky='we', padx=5, pady=3)
        text.grid(row=5, column=1, sticky='we', padx=5, pady=3)
        return parent

    def _fillInUI_Prody(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

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

        #
        # Copyright
        # 
        row += 1
        self.ui_banner = self._team_logo(self.canvas)
        self.ui_banner.grid(row=row, columnspan=2, sticky='news', padx=5, pady=5)

    def _fillInUI_Gaussian(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

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

        #
        # Copyright
        # 
        row += 1
        self.ui_banner = self._team_logo(self.canvas)
        self.ui_banner.grid(row=row, columnspan=2, sticky='news', padx=5, pady=5)

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

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        chimera.extension.manager.deregisterInstance(self)
        self.destroy()
    
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


class NormalModesResultsDialog(ModelessDialog):

    buttons = ('Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, parent=None, controller=None, *args, **kwarg):
        # GUI init
        self.title = 'Normal Modes Results'
        self.parent = parent
        self.controller = controller
        # Fire up
        ModelessDialog.__init__(self, resizable=False)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def fillInUI(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        chimera.extension.manager.deregisterInstance(self)
        self.destroy()

    def populate_data(self, frequencies=None):
        pass

    def plot_vectors(self, vectors=None):
        pass


class NormalModesMovieDialog(ModelessDialog):

    buttons = ('Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, parent=None, controller=None, *args, **kwarg):
        # GUI init
        self.title = 'Normal Modes Results'
        self.parent = parent
        self.controller = controller
        # Fire up
        ModelessDialog.__init__(self, resizable=False)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def fillInUI(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        # self.destroy()


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