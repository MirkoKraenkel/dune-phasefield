#! /usr/bin/env python
import subprocess,sys,Tkinter
sys.path.append('../../scripts/')
from Tkinter import *
from pfmaker import *

root=Tk()
app=PhasefieldMaker( root )
root.mainloop()



