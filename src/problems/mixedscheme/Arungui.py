#! /usr/bin/env python
import subprocess,sys,Tkinter
sys.path.append('../../scripts/')
from Tkinter import *
from pfArunner import *

root=Tk()
app=PhasefieldARunner( root )
root.mainloop()



