TOP	= ../..
include $(TOP)/mk/common.make

MODULE	= nmod
MODDIR	= $(datadir)/$(MODULE)
SRCS	= ChimeraExtension.py MMTK2Molecule.py NormalModesTable.py \
		OnPicker.py __init__.py base.py cmdline.py \
		confChange.py gui.py nmodMMTKinter.py plotdialog.py
OBJS	= $(SRCS:.py=.pyc)
CLEAN	= $(OBJS)

all: $(OBJS)

install: all
	-mkdir -p $(MODDIR)
	$(RSYNC) $(SRCS) $(OBJS) $(DATA) $(MODDIR)
