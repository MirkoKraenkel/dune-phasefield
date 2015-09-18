#!/usr/bin/env python
import subprocess, sys,argparse

def call( arg ):
    subprocess.call( arg ,shell=True)

parser = argparse.ArgumentParser()
parser.add_argument('executable')
parser.add_argument('dirname')
args=parser.parse_args()

topdir=args.dirname+'/'
subdir=topdir+args.executable
#subprocess.call('cp '+args.executable+' foo', shell=True)
call('mkdir '+topdir )
call('mkdir '+subdir)
call('cp '+args.executable+' '+subdir)
call('cp -r paramFiles/ '+subdir)
call('cp -r ../parameterFiles/ '+topdir)
call('cp -r ../macrogrids/ '+topdir)
call('cp -r ../InputFiles/ '+topdir) 


