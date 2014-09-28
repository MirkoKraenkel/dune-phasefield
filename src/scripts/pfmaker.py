import subprocess,sys,pickle
from Tkinter import *
terminal="xfce4-terminal"

#==================================
#the make class
#==================================
class PhasefieldMaker:
  def __init__(self, master):
    entryframe=Frame(master)
    entryframe.pack()
    listframe=Frame(master)
    listframe.pack()
    messageframe=Frame(master)
    messageframe.pack(side=BOTTOM)
    self.cxxflags={ 'debug':'-Wfatal-errors -g -Wall' , 'full' :''}
    self.compiled=pickle.load(open("compiled.p","rb"))
    self.programms=[]
    progfile=open('progs.txt')
    for line in progfile:
      self.programms.append(line.rstrip('\n'))
    self.programms
    progfile.close()
    self.listbox = Listbox( listframe , selectmode=MULTIPLE)
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
    self.msgvar.set(argstring)
  def composeString( self ):
    flag=''
    if self.var.get == 'debug':
      flag='CXXFLAGS="'+ self.cxxflags[self.var.get()]+'"'
    polOrder=int( self.polord.get())
    idxs = self.listbox.curselection()
    for ind in idxs:
      index=int(ind)
      p=self.programms[index]
      outfile=p+'_make.out'
      execstring= 'make '+flag+' POLORDER='+ str(polOrder)+' '+p+' &>'+outfile
      self.updateMsg(execstring)
      subprocess.call(['make clean'],shell=True)
      compiled=self.makecall(( execstring , outfile ))
      if compiled == 0:
        realname=p+'_'+str( polOrder )
        subprocess.call(['mv '+p+' '+realname],shell=True)
        self.compiled.append( realname ) 
        self.updateMsg('succeeded')
    pickle.dump( self.compiled, open("compiled.p","wb"))
  def makecall(self, stringtuple):
    execstring=stringtuple[0]
    print(execstring)
    outfile=stringtuple[1]
    #subprocess.call([terminal+' -e "'+execstring+'"' ],shell=True)  
    self.updateMsg( execstring )
    c = subprocess.call([execstring], shell=True)
    if c == 0:
      subprocess.call(['rm '+outfile],shell=True)
    else:
      res=self.scanoutfile(outfile)
      subprocess.call([terminal+' -e "vi '+outfile+'"' ],shell=True)  
      subprocess.call([terminal+' -e "vi +'+res[1]+' '+res[0]+'"' ],shell=True) 
    return c
  def makeit(self):
    self.composeString()
  def clean(self):
    self.updateMsg('clean!')
    subprocess.call(['make clean'],shell=True)
    todelete=pickle.load(open("compiled.p","rb"))
    for p in todelete:
      print(p)
      subprocess.call(['rm '+p],shell=True)
    self.compiled=[]
    pickle.dump( self.compiled, open("compiled.p","wb"))
  def scanoutfile(self,filename):
    myfile=open(filename)
    for line in myfile:
      if line.find('error') !=-1:
        splitline=line.split(":")
        break
    return (splitline[0], splitline[1])
