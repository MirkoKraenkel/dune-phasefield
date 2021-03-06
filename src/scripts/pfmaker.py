import subprocess,sys,pickle, os.path
from Tkinter import *
terminal="xfce4-terminal"

#==================================
#the make class
#==================================
class PhasefieldMaker:
  
  def __init__(self, master):
    flagframe=Frame(master)
    flagframe.pack()
    listframe=Frame(master)
    listframe.pack()
    listframe2=Frame(master)
    listframe2.pack()
    entryframe=Frame(master)
    entryframe.pack()
    messageframe=Frame(master)
    messageframe.pack(side=BOTTOM)
    self.cxxflags={ 'debug':'-Wfatal-errors -g -Wall' , 'full' :''}
    self.flags=['POLORDER','GRIDDIM','GRIDTYPE','PROBLEM']
    self.defaults=['1','2','ALUGRID_CUBE','2']
    MYPATH = './flags.p'
    if os.path.isfile(MYPATH):
      self.vals=pickle.load(open("flags.p","rb"))
    else:
      self.vals=dict(zip(self.flags,self.defaults))
    COMPILEDPATH='./compiled.p'
    if os.path.isfile(COMPILEDPATH):
      self.compiled=pickle.load(open("compiled.p","rb"))
    else:
      self.compiled=dict()
    self.programms=[]
    progfile=open('progs.txt')
    for line in progfile:
      self.programms.append(line.rstrip('\n'))
    self.programms
    progfile.close()
    self.listbox = Listbox( listframe , selectmode=MULTIPLE, exportselection=0)
    self.listbox.pack(side=LEFT)
    for item in self.programms:
      self.listbox.insert(END, item)
    self.listbox2 = Listbox( listframe2 , selectmode=MULTIPLE, exportselection=0)
    self.listbox2.pack(side=LEFT)
    self.thermlist=[]
    if os.path.isfile('./thermo.txt'):
      thermfile=open('thermo.txt')
      for line in thermfile:
        self.thermlist.append(line.rstrip('\n'))
      thermfile.close()
    for item in self.thermlist:
        self.listbox2.insert(END, item)
    self.var = StringVar() 
    for item in self.cxxflags:
      var=IntVar()
      chck =Radiobutton( entryframe, text=item ,variable=self.var, value=item).pack(side=LEFT)
    self.make=Button( entryframe, text="Make!", command=self.makeit)
    self.make.pack(side=RIGHT)
    self.flagentries={}
    self.makeform(flagframe,self.vals,self.flagentries) 
    self.msgvar=StringVar()
    self.msgvar.set('welcome')
    mymsg = Message( master, textvariable=self.msgvar,width=250).pack()
    self.cleanbutton= Button( master, text="Clean!", command=self.clean)
    self.cleanbutton.pack(side=BOTTOM)
  
  
  #create form from flags
  def makeform( self , root ,fields,entries):
    ##entries={}
    for field in fields:
      row=Frame(root)
      lab=Label(row,text=field)
      ent=Entry(row)
      row.pack(side=TOP,fill=X)
      lab.pack(side=LEFT)
      ent.pack(side=RIGHT,fill=X,expand=YES)
      ent.insert(0,fields[field])
      entries.update({field : ent})
    #return entries
  
  
  def updateMsg( self, argstring):
    self.msgvar.set(argstring)
  
  
  #create string with compile flags an pass it to make (method name makes no sense)
  def composeString( self ):
    flag=''
    
    if self.var.get == 'debug':
      flag='CXXFLAGS="'+ self.cxxflags[self.var.get()]+'"'
    
    flag2=''
    flag3=''
    addstring=''
    outdict={} 
    for key in self.flagentries:
      flag2+=key+'='+self.flagentries[key].get()+' '
      outdict.update({key:self.flagentries[key].get()})
    
    flag+=flag2
    idxs = self.listbox.curselection()
    therms = self.listbox2.curselection() 
    
    for ind in idxs:
      for term in therms:
        index=int(ind)
        index2=int(term)
        thermflag='THERMO='+str(index2+1)
        p=self.programms[index]
        outfile=p+'_make.out'
        execstring= 'make '+flag+' '+thermflag+' '+p + ' &>'+outfile
        self.updateMsg(execstring)
        subprocess.call(['make clean'],shell=True)
        compiled=self.makecall(( execstring , outfile ))
      
        if compiled == 0:
            realname=p+'_'+self.thermlist[index2]+'_P'+str(self.flagentries['POLORDER'].get())+'_PROB'+str(self.flagentries['PROBLEM'].get())
            subprocess.call(['mv '+p+' '+realname],shell=True)
            self.compiled.append( realname ) 
            self.updateMsg('succeeded')
    pickle.dump( outdict, open("flags.p", "wb")) 
    pickle.dump( self.compiled, open("compiled.p","wb"))
  
  
  #call make in and scan outputfiles, in case of errors open output- and sourcefiles
  def makecall(self, stringtuple):
    execstring=stringtuple[0]
    print(execstring)
    outfile=stringtuple[1]
    self.updateMsg( execstring )
    c = subprocess.call([execstring], shell=True)
    if c == 0:
      subprocess.call(['rm '+outfile],shell=True)
    else:
      res=self.scanoutfile(outfile)
      subprocess.call([terminal+' -e "vi '+outfile+'"' ],shell=True)  
      subprocess.call([terminal+' -e "vi +'+res[1]+' '+res[0]+'"' ],shell=True) 
    return c
  
  #forward call(redesign...) 
  def makeit(self):
    self.composeString()
  
  #clean
  def clean(self):
    self.updateMsg('clean!')
    subprocess.call(['make clean'],shell=True)
    todelete=pickle.load(open("compiled.p","rb"))
    for p in todelete:
      print(p)
      subprocess.call(['rm '+p],shell=True)
    self.compiled=[]
    pickle.dump( self.compiled, open("compiled.p","wb"))
  
  #find error messages in compiler output
  def scanoutfile(self,filename):
    myfile=open(filename)
    splitline=[]
    for line in myfile:
      if line.find('error') !=-1:
        splitline+=line.split(":")
        break
    print(splitline)
    return (splitline[0], splitline[1])
