import chimera
import Tkinter as tk

from chimera.baseDialog import ModelessDialog

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
        pass

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

    def fillInUI_prody(self,parent):
        self.parent = parent
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')


