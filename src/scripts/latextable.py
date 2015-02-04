import os, subprocess, sys
infile = (sys.argv[1])
outfile = (sys.argv[2])
fin=open(infile)
fout=open(outfile, 'w')
fout.write("\\begin{tabular}{|")
tabheader=''

d=fin.readlines()
indices=[1,4,5,6,8]
for i in indices:
  tabheader+="c|"
tabheader+="}\n"
fout.write(tabheader)
fout.write("\\cline\n")
myheader=''
for i in indices:
  col = d[0].split()[1+i]
  myheader+=col+" & "
myheader=myheader[:-2]
myheader+="\n"
fout.write(myheader)
fout.write("\\cline\n")
for row in d[1:]:
  mystring=''
  for i in indices:
    col = row.split()[i]
    mystring+="$"+col+"$"+" & "
  mystring=mystring[:-2]
  mystring+="\n"
  fout.write(mystring)
fout.write("\\end{tabular}")
