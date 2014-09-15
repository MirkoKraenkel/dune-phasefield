#! /usr/bin/env python
import subprocess,sys
from Tkinter import *
terminal="xfce4-terminal"

#==================================
#the make class
#==================================
class PhasefieldMaker:
  def __init__(self, master):
    listframe=Frame(master)
    listframe.pack()
    entryframe=Frame(master)
    entryframe.pack()
    messageframe=Frame(master)
    messageframe.pack(side=BOTTOM)
    self.cxxflags={ 'debug':'-Wfatal-errors -g -Wall' , 'full' :''}
    self.programms=[ 'phasefield', 'phasefieldcoupling','phasefieldfd', 'phasefieldfdrho', 'phasefieldfdlambda' , 'phasefieldmfree', 'phasefieldoem']
    self.listbox = Listbox( listframe )
    self.listbox.pack(side=LEFT)
    for item in self.programms:
      self.listbox.insert(END, item)
    self.var = StringVar() 
    for item in self.cxxflags:
      var=IntVar()
      chck =Radiobutton( entryframe, text=item ,variable=self.var, value=item).pack(side=LEFT)
    self.make=Button( entryframe, text="Make!", command=self.makeit)
    self.make.pack(side=RIGHT)
    self.polord = Entry( entryframe)
    self.polord.insert(0,'PolOrder')
    self.polord.pack(side=RIGHT )
    self.msgvar=StringVar()
    self.msgvar.set('welcome')
    self.mymsg = Message( master, textvariable=self.msgvar,width=250).pack()
    self.cleanbutton= Button( master, text="Clean!", command=self.clean)
    self.cleanbutton.pack(side=BOTTOM)
  def updateMsg( self, argstring):
    print("upadate")
    self.msgvar.set(argstring)
  def composeString( self ):
    flag = self.cxxflags[self.var.get()]
    idxs = self.listbox.curselection()
    index=int(idxs[0])
    polOrder=int( self.polord.get())
    p=self.programms[index]
    outfile=p+'_make.out'
    execstring= 'make CXXFLAGS="'+flag+'" POLORDER='+ str(polOrder)+' '+p+' &>'+outfile
    self.updateMsg(execstring)
    return ( execstring , outfile )
  def makecall(self, stringtuple):
    execstring=stringtuple[0]
    outfile=stringtuple[1]
    c = subprocess.call([execstring], shell=True)
    print( c )
    if c == 0:
      subprocess.call(['rm '+outfile],shell=True)
    else:
      subprocess.call([terminal+' -e "vi '+outfile+'"' ],shell=True) 
  def makeit(self):
    self.makecall(self.composeString())
  def clean(self):
    self.updateMsg('clean!')
    subprocess.call(['make clean'],shell=True)
root=Tk()
app=PhasefieldMaker( root )

root.mainloop()



