#!/usr/bin/env python

# Create objects to manipulate - they hold your data
#
import sys
from ImageD11 import peakmerge, indexing, transformer
mypeakmerger = peakmerge.peakmerger()
mytransformer = transformer.transformer()
myindexer = indexing.indexer()
#
# Your work starts here:
#
filename = sys.argv[1]
outfile =  sys.argv[2]
parfile =  sys.argv[3]
mytransformer.loadfiltered(  filename )
mytransformer.loadfileparameters( parfile )
mytransformer.compute_tth_eta( )
mytransformer.getcolumn(  'tth' )
mytransformer.getcolumn(  'eta' )
mytransformer.addcellpeaks( )
mytransformer.computegv( )
mytransformer.savegv(  outfile )
