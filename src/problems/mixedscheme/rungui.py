#! /usr/bin/env python
import subprocess,sys,Tkinter
sys.path.append('../../scripts/')
from Tkinter import *
from pfrunner import *

root=Tk()
app=PhasefieldRunner( root )
root.mainloop()



