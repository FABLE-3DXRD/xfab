#!/usr/bin/env python

'''
grainspotter_loop

Runs GrainSpotter until no more good seed grains are found
Should enter one variable: the full path to the script's INI file
Like: /home/jules/sp8_2013/inis/script.ini
Must be run in folder where n00.flt is located
'''

from __future__ import absolute_import
from __future__ import print_function
import sys
import os
import shutil
import re
from string import split,digits,letters#,replace


GS='grainspotter'
GFF2UBI='gff_to_ubi.py'
RUNMAKEMAP='makemap_all.py'
MAKEGVE='make_gve.py'



# read and check input file
if len(sys.argv) == 1:
    print('Provide path to SCRIPT.ini')
elif len(sys.argv) != 2:
    print('Wrong number of arguments entered')
    sys.exit()
else:
    pathtoINI = sys.argv[1]
    print('\nCorrect number of arguments entered -- using %s \n' %pathtoINI)
    
f = open(pathtoINI,'r')
input = f.readlines()
f.close()

INIFILE=False
for line in input:
    if 'INIFILE' in line:
        INIFILE = eval(split(line,'=')[1])
        break
if INIFILE==False:
    print('Missing INIFILE variable - update ini file and restart')
else:
    print(' GrainSpotter ini file: %s' %INIFILE)
    f = open(INIFILE,'r')
    input = f.readlines()
    f.close()
    for line in input:
        if 'omegarange' in line:
            minome = eval(split(line)[1])
            maxome = eval(split(line)[2])
            break
    if minome < -180 or maxome > 180:
        print('\n OUT OF RANGE: omegarange (%s) must be within -180 to +180.' % INIFILE)
        sys.exit()
    print(' Omega range: %0.2f to %0.2f' %(minome,maxome))

f = open(pathtoINI,'r')
input = f.readlines()
f.close()
    
fltfile=False
for line in input:
    if 'fltfile' in line:
        fltfile = eval(split(line,'=')[1])
        break
if fltfile==False:
    print('Missing fltfile variable - update ini file and restart')
else:
    print(' GrainSpotter input peaks: %s' %fltfile)
    
ubifile=False
for line in input:
    if 'ubifile' in line:
        ubifile = eval(split(line,'=')[1])
        prefix  = split(ubifile,'.')[0].translate(None,digits)
        itemnum = eval(split(ubifile,'.')[0].translate(None,letters))
        break
if ubifile==False:
    print('Missing ubifile variable - update ini file and restart')
else:
    print(' GrainSpotter output grains: %s' %ubifile)
    
NUMFOUND=False
for line in input:
    if 'NUMFOUND' in line:
        NUMFOUND = eval(split(line,'=')[1])
        break
if NUMFOUND==False:
    print('Missing NUMFOUND variable - update ini file and restart')
else:
    print(' Continue running GrainSpotter while the number of new grains is at least: %s' %NUMFOUND)
    
PARS=False
for line in input:
    if 'PARS' in line:
        PARS = eval(split(line,'=')[1])
        break
if PARS==False:
    print('Missing PARS variable - update ini file and restart')
else:
    print(' Full path to par file is: %s' %PARS)
    
ICUTS=False
for line in input:
    if 'ICUTS' in line:
        ICUTS = eval(split(line,'=')[1])
        break
if ICUTS==False:
    print('Missing ICUTS variable - update ini file and restart')
else:
    print(' Makemap initial cuts: %s' %ICUTS)
    
FCUTS=False
for line in input:
    if 'FCUTS' in line:
        FCUTS = eval(split(line,'=')[1])
        break
if FCUTS==False:
    print('Missing FCUTS variable - update ini file and restart')
else:
    print(' Makemap final cuts: %s' %FCUTS)
    
SLIMIT=False
for line in input:
    if 'SLIMIT' in line:
        SLIMIT = eval(split(line,'=')[1])
        break
if SLIMIT==False:
    print('Missing SLIMIT variable - update ini file and restart')
else:
    print(' Makemap stop limit, number of new reflections per grain: %s' %SLIMIT)

MINPEAKS=False
for line in input:
    if 'MINPEAKS' in line:
        MINPEAKS = eval(split(line,'=')[1])
        break
if MINPEAKS==False:
    print('Missing MINPEAKS variable - update ini file and restart')
else:
    print(' Makemap minimum number of peaks: %s' %MINPEAKS)

CRYSYMM=False
for line in input:
    if 'CRYSYMM' in line:
        CRYSYMM = eval(split(line,'=')[1])
        break
if CRYSYMM==False:
    print('Missing CRYSYMM variable - update ini file and restart')
else:
    print(' The crystal symmetry is: %s' %CRYSYMM)    
    
EULERSTEP=False
for line in input:
    if 'EULERSTEP' in line:
        ES = split(line,'=')[1]
        EULERSTEP = eval(ES.replace('(','[').replace(')',']').replace(' ',','))
        break
if EULERSTEP==False:
    print('Missing EULERSTEP variable - update ini file and restart')

NSIGMAS=False
for line in input:
    if 'NSIGMAS' in line:
        NS = split(line,'=')[1]
        NSIGMAS = eval(NS.replace('(','[').replace(')',']').replace(' ',','))
        break
if NSIGMAS==False:
    print('Missing NSIGMAS variable - update ini file and restart')

    
# Setup for loop
testCon=0
testLim=3
loopCounter=0
numberofgrainsbefore=0
numberofgrains=-1
spottercounter=-1


while testCon < testLim:
    grains_in = '%s%0.2d.ubi' %(prefix,itemnum+loopCounter)
    grains_out = '%s%0.2d.ubi' %(prefix,itemnum+loopCounter)
    pks_in = '%s%0.2d.flt' %(prefix,itemnum-1+loopCounter)    
    pks_out = '%s%0.2d.flt' %(prefix,itemnum+loopCounter)
    log_out = '%s%0.2d.log' %(prefix,itemnum+loopCounter)
    gff_out = '%s%0.2d.gff' %(prefix,itemnum+loopCounter)
    
    print('\nConvert peaks to g-vectors')
    command = '%s %s rest.gve %s' %(MAKEGVE,pks_in,PARS)
    #print command
    os.system(command)

    if numberofgrains < NUMFOUND:
        spottercounter = spottercounter + 1        
    if spottercounter >= len(EULERSTEP) or spottercounter >= len(NSIGMAS):
        sys.exit(0)

    eulerstep = EULERSTEP[spottercounter]
    nsigmas = NSIGMAS[spottercounter]
    print('\nRunning GrainSpotter')
    print('spottercounter is: %i' %spottercounter)
    print('Now using eulerstep: %0.1f and nsigmas %0.1f' %(eulerstep,nsigmas))
    
    ## REPLACES cat COMMAND IN ORDER TO FUNCTION ON WINDOWS
    ## mkak 14.03.2016
    NEWINIFILE = 'grainspotter_current.ini'
    shutil.copyfile(INIFILE, NEWINIFILE)
    with open(NEWINIFILE, 'r') as sources:
        lines = sources.readlines()
    with open(NEWINIFILE, 'w') as sources:
        for line in lines:
            sources.write(re.sub(r'_EULERSTEP_', '%0.6f' % eulerstep, line))
    with open(NEWINIFILE, 'r') as sources:
        lines = sources.readlines()
    with open(NEWINIFILE, 'w') as sources:
        for line in lines:
            sources.write(re.sub(r'_NSIGMAS_', '%0.6f' % nsigmas, line))

    
    command = '%s grainspotter_current.ini > spotterOutTemp.txt' %GS
    os.system(command) 
    f = open('spotter_out.log','r')
    lines = f.readlines()
    f.close()
    print(lines[0], 'GrainSpotter done')
    numberofgrainsbefore = numberofgrains
    numberofgrains = eval(lines[0].translate(None,letters))
    slimit = SLIMIT*numberofgrains    
    
    
    ## REPLACES os.rename WITH shututil.move IN ORDER TO FUNCTION THE SAME ON ALL OS
    ## mkak 17.03.2016
    shutil.move("spotter_out.gff",gff_out)
    shutil.move("spotter_out.log",log_out)
    os.remove('spotter_out.ubi')

    if numberofgrains != 0:
        command = '%s %s %s' %(GFF2UBI,gff_out,grains_in)
        os.system(command)
        print('\nRunning makemap.py')
        command = '%s %s %s %s %f %f %i %i %s %s' %(RUNMAKEMAP,grains_in,pks_in,pks_out,ICUTS,FCUTS,slimit,MINPEAKS,PARS,CRYSYMM)
        os.system(command)
        
        if os.path.exists(grains_in):
            testCon = 0
            loopCounter = loopCounter + 1
        else:
            testCon = testCon + 1
    else:
        print('No grains found')
        os.remove(gff_out)
        os.remove(log_out)
        testCon = testCon + 1
        
os.remove('spotterOutTemp.txt')
print('\nOK\n')

        