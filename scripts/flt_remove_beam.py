#!/usr/bin/env python



##########################################################################
# IMPORT PYTHON PACKAGES
from __future__ import absolute_import
from __future__ import print_function
import os
import argparse
import textwrap
from six.moves import range

##########################################################################
def sortFLT(FLTfile,newFLTfile,args):

    f=open(FLTfile,'r')
    fw=open(newFLTfile,'w')

    for i in range(2):
        line=f.readline()
        fw.write(line)
    titles = line.split()
   
    eof = 0
    while not(eof):
        line = f.readline()
        s = line.split()
        if len(s)>0:
       
            if float(s[1]) < args.bmaxx and float(s[1]) > args.bminx and float(s[0]) < args.bmaxy and float(s[0]) > args.bminy:
                a = 11
            else:
                fw.write('%s' % line)

        else:
            eof = 1

    f.close()
    fw.close()
    
    print('\nCreated file %s.' % (newFLTfile))


    return()

##########################################################################
def main(args):

    if not os.path.exists(args.fltfile):
    	raise IOError('File %s does not exist.' % args.fltfile)

    newfltfile =  '%s_%s.flt' % (args.fltfile.split('.')[0],'noBEAM')
    if os.path.exists(newfltfile) and not args.force:
    	raise IOError('File %s already exists.' % newfltfile)


    if args.bmaxx < args.bminx:
        newmax = args.bminx
        args.bminx = args.bmaxx
        args.bmaxx = newmax
        if args.verbose:
            print('Redefined limits of beamspot: x min and max were switched.')
    if args.bmaxy < args.bminy:
        newmax = args.bminy
        args.bminy = args.bmaxy
        args.bmaxy = newmax
        if args.verbose:
            print('Redefined limits of beamspot: y min and max were switched.')

    sortFLT(args.fltfile,newfltfile,args)

    return()

##########################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='flt_remove_beam.py',
        description='Removes peaks in within the near-field beam spot from flt file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
             examples:
 
             $ 
             $ flt_remove_beam.py peaks_t100.flt 250 1900 935 1125
             ''')
         )

    parser.add_argument(
        'fltfile', type=str, help='flt file'
        )
    parser.add_argument(
        'bminx', type=int, help='beam spot limits: min. x'
        )        
    parser.add_argument(
        'bmaxx', type=int, help='beam spot limits: max. x'
        )        
    parser.add_argument(
        'bminy', type=int, help='beam spot limits: min. y'
        )        
    parser.add_argument(
        'bmaxy', type=int, help='beam spot limits: max. y'
        )        
    parser.add_argument(
        '-f', '--force', default=False, action='store_true',
        help='force removal of report file for overwriting'
        )  
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='report progress in terminal'
        )
    args = parser.parse_args()
    main(args)
