#!/opt/paraviewopenfoam410/bin/pvpython

from paraview.simple import *
import sys
import os


if len(sys.argv) > 1:
    inputFilename = os.getcwd() + '/' +  sys.argv[1]
    outputFilename = os.getcwd() + '/' + sys.argv[1].split(".")[0] + ".stl"

    pvFilename = ''.join(outputFilename.split(".")[:-1]) + "0.stl"

    #read VRML
    reader = VRMLReader(FileName=inputFilename)

    #write STL
    writer = PSTLWriter(reader, FileName=outputFilename)
    writer.UpdatePipeline()
    os.rename(pvFilename, outputFilename)

else:
    print 'Usage: vrml2stl <input.wrl>'
