#! /usr/bin/env python
import subprocess,sys
from Tkinter import *
terminal="xfce4-terminal"

#==================================
#the make class
#==================================
class PhasefieldRunner:
  def __init__(self, master):
    listframe=Frame(master)
    listframe.pack()
    paramframe=Frame(master)
    paramframe.pack()
    entryframe=Frame(master)
    entryframe.pack()
    messageframe=Frame(master)
    messageframe.pack(side=BOTTOM)
    self.programms=[ 'phasefield', 'phasefieldcoupling','phasefieldfd', 'phasefieldfdrho', 'phasefieldfdlambda' , 'phasefieldmfree', 'phasefieldoem']
    self.params=['delta','mu1','mu2','Endtime', 'deltaT','out path','startLevel','runs']
    self.defaults=['0.01','0.1','0.1','1','1e-3','./data/','2', '1']
    self.paramEntries=self.makeform(paramframe,self.params)
    self.listbox = Listbox( listframe )
    self.listbox.pack(side=LEFT)
    for item in self.programms:
      self.listbox.insert(END, item)
    self.make=Button( entryframe, text="Run", command=self.runit)
    self.make.pack(side=RIGHT)
  def makeform( self , root ,fields):
    entries={}
    count=0
    for field in fields:
      row=Frame(root)
      lab=Label(row,text=field)
      ent=Entry(row)
      row.pack(side=TOP,fill=X)
      lab.pack(side=LEFT)
      ent.pack(side=RIGHT,fill=X,expand=YES)
      ent.insert(0,self.defaults[count])
      entries.update({field : ent.get()})
      count+=1
    return entries
  def composeString( self ):
    idxs = self.listbox.curselection()
    index=int(idxs[0])
    p=self.programms[index]
    outfile=self.paramEntries['out path']
    execstring ='./'+p+' parameter_gui delta:'+self.paramEntries['delta']\
    +' mu1:'+self.paramEntries['mu1']+' mu2:'+self.paramEntries['mu2'] \
    +' endTime:'+self.paramEntries['Endtime'] \
    +' deltaT:'+self.paramEntries['deltaT']\
    +' startLevel:'+self.paramEntries['startLevel']\
    +' path:'+outfile+'_'+p
    print( outfile ) 
    return ( execstring , outfile )
  def runcall(self, stringtuple):
    execstring=stringtuple[0]
    outfile=stringtuple[1]
    c = subprocess.call([execstring], shell=True)
  def runit(self):
    self.runcall(self.composeString())

root=Tk()
app=PhasefieldRunner( root )

root.mainloop()



