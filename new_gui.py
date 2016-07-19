import chimera
import Tkinter as tk

from chimera.baseDialog import ModelessDialog

from core import Controller

class NormalModesExtension(ModelessDialog):
    buttons = ('OK','Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, *args, **kwarg):
        # GUI init
        self.title = 'Plume Normal Modes'
        # Fire up
        ModelessDialog.__init__(self,reseizable=False)
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

        self.input_frame = tk.LabelFrame(self.canvas, text='Select mode', padx=5, pady=5)
        self.input_frame.pack(expand=True, fill='x')
        self.input_choice_frame = tk.Frame(self.input_frame)
        self.input_choice_frame.grid(row=0)
        self.input_choice = tk.StringVar()
        self.input_choice.set('prody')
        self.input_choice_prody = tk.Radiobutton(self.input_choice_frame, variable=self.input_choice,
                                                     text='ProDy', value='prody')
        self.input_choice_gaussian = tk.Radiobutton(self.input_choice_frame, variable=self.input_choice,
                                                     text='Gaussian', value='gaussian')
        self.input_choice_prody.pack(side='left')
        self.input_choice_gaussian.pack(side='left')
        self.input_choice_prody.select()

    def Apply(self):
        """
        Default! Triggered action if you click on an Apply button
        """
        if self.modes_dialog is None:
            self.modes_dialog = NormalModesConfigDialog(self)
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
        # self.destroy()


class NormalModesConfigDialog(ModelessDialog):

    buttons = ('OK','Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, parent=None, *args, **kwarg):
        # GUI init
        self.title = 'Calc Normal Modes'
        self.parent = parent
        slef.load_method()
        # Fire up
        ModelessDialog.__init__(self,reseizable=False)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def fillInUI_Prody(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

    def fillInUI_Gaussian(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

    def Apply(self):
        """
        Default! Triggered action if you click on an Apply button
        Change in core for apply_prody or apply_gaussian
        """
        self.controller = Controller(self)
        self.vibrations = self.controller.run()
        if self.results_dialog is None:
            self.results_dialog = NormalModesResultsDialog(self)
        self.results_dialog.enter()
        self.results_dialog.fillInData()
        if self.movie_dialog is None:
            slef.movie_dialog = NormalModesMovieDialog(self)
        self.results_dialog.enter()

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
        # self.destroy()

    def load_method(self):
        if self.parent.input_choice == 'prody'
            slef.fillInUI = self.fillInUI_Prody
        else:
            self.fillInUI = self.fillInUI_Gaussian

class NormalModesResultsDialog(ModelessDialog):

    buttons = ('Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, parent=None, *args, **kwarg):
        # GUI init
        self.title = 'Normal Modes Results'
        self.parent = parent
        # Fire up
        ModelessDialog.__init__(self,reseizable=False)
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

    def fillInData(self):
        pass

class NormalModesMovieDialog(ModelessDialog):

    buttons = ('Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, parent=None, *args, **kwarg):
        # GUI init
        self.title = 'Normal Modes Results'
        self.parent = parent
        # Fire up
        ModelessDialog.__init__(self,reseizable=False)
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