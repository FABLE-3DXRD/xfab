#!/usr/bin/env python
from __future__ import print_function

import glob, os

def lintit(infile,outfile):
    print(infile,outfile)
    os.system("pylint %s > %s"%(infile,outfile))


pyfiles = glob.glob(os.path.join("..","xfab","*.py"))

for f in pyfiles:
    outf = os.path.split(f)[-1] + ".lint"
    if not os.path.exists(outf) :
        lintit(f,outf)
        continue
    if os.stat(f).st_mtime > os.stat(outf).st_mtime:
        lintit(f,outf)
        continue

