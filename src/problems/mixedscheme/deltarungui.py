#! /usr/bin/env python
import subprocess,sys,Tkinter
sys.path.append('../../scripts/')
from Tkinter import *
from pfdeltarunner import *

root=Tk()
app=PhasefieldDeltaRunner( root )
root.mainloop()



