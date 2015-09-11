#! /usr/bin/env python
import subprocess,sys,pickle,time
from Tkinter import *
terminal="xfce4-terminal"

#==================================
#the Runner class
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
    self.params=[ 'paramfile','phasefield.delta','phasefield.mu1','phasefield.mu2','phasefield.endTime', 'fixedTimeStep',
                  'fem.prefix','phasefield.startLevel','phasefield.reactionrate','phasefield.alpha',
                  'phasefield.beta','phasefield.A', 'phasefield.rhofactor','phasefield.addvisc',
                  'phasefield.muvisc','phasefield.acpenalty','phasefield.gravity','phasefield.timestepfactor','phasefield.rescale']
    self.defaults=['*.log', '0.01','0.1','0.1','1','0','','0','1','0.','1','0.002','0',1,0., 1,0.,1,'false']
    
    self.paramEntries=self.makeform(paramframe,self.params)
    
    self.listbox = Listbox( listframe, width=30 )
    self.listbox.pack(side=LEFT)
    
    self.programms=pickle.load(open("compiled.p","rb"))
    
    for item in self.programms:
      self.listbox.insert(END, item)
    
    self.make=Button( entryframe, text="Run", command=self.runit)
    self.make.pack(side=RIGHT)
    self.scrvar=IntVar()
    self.rdchck=IntVar()
    readcheck=Radiobutton(entryframe,text='read from file', variable=self.rdchck,value=1).pack()
    readcheck=Radiobutton(entryframe,text='read from form', variable=self.rdchck,value=0).pack()

    scr=Radiobutton(entryframe,text='screen',variable=self.scrvar,value=1).pack() 
    scr=Radiobutton(entryframe,text='show',variable=self.scrvar,value=0).pack() 
  
###create form from params#####################
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
      entries.update({field : ent})
      count+=1
    
    return entries
##############################################

  def composeFileString( self ):

    myday=time.strftime("%d_%m")
    idxs = self.listbox.curselection()
    index=int(idxs[0])
    p=self.programms[index]
    outfile='/disk1/Data_'+myday+'/'+self.paramEntries['fem.prefix'].get()
    paramstring=''
    paramstring+='fem.prefix:'+outfile+'_'+p+' '
    execstring ='./'+p+' '+paramstring+' paramfile:'+self.paramEntries['paramfile'].get()
    return ( execstring , outfile )
  
  #compose string passed to the command line
  def composeFormString( self ):
    
    myday=time.strftime("%d_%m") 
    idxs = self.listbox.curselection()
    index=int(idxs[0])
    p=self.programms[index]
    outfile='/disk1/Data_'+myday+'/'+self.paramEntries['fem.prefix'].get()
    paramstring=''
    
    for key in self.paramEntries:
      
      if key =='fem.prefix':
        
        paramstring+=key+':'+outfile+'_'+p+' '
      
      else:
       
       paramstring+=key+':'+self.paramEntries[key].get()+' '
    
    #execstring ='/opt/openmpi/bin/mpiexec -n 4 '+p+' '+paramstring+'paramfile:parameter_gui'
    #execstring ='/opt/openmpi/bin/mpiexec -n 4  xterm -e gdb --args '+p+' '+paramstring+'paramfile:parameter_gui'
    execstring ='./'+p+' '+paramstring+'paramfile:../parameterFiles/parameter'
    #execstring ='./'+p+' '+paramstring+'paramfile:../parameterFiles/parameter &> '+outfile+p+ '.out'
    #

    #execstring ='gdb --args '+p+' '+paramstring+'paramfile:parameter_gui '
    return ( execstring , outfile )
  
  #call the programm with parameter string directly or in a screen session
  def runcall(self, stringtuple):
    
    execstring=stringtuple[0]
    print(execstring)
    outfile=stringtuple[1]
    print( outfile )
    print( self.scrvar.get()) 
    
    if self.scrvar.get() == 0: 
      
      c = subprocess.call([execstring], shell=True)
    
    else:
      
      c = subprocess.call([terminal+' -e "screen '+execstring+'"'], shell=True)
      #subprocess.call(['screen -d'],shell=True)
  
  #run the programm
  def runit(self):

    if self.rdchck.get()==1:

        self.runcall(self.composeFileString())

    else:

        self.runcall(self.composeFormString())
