import os
import sys
import argparse

sys.path.insert(0, "/home/redo/Sources/GridDataFormats/")

from gridData import Grid

import numpy as np

import basicutils

if __name__ == "__main__":

     parser = argparse.ArgumentParser()
     parser.add_argument("-f1","--inputfile1", help="Specify DX input file", required=True, 
             type=str, default="")
     parser.add_argument("-f2","--inputfile2", help="Specify DX input file", required=True, 
             type=str, default="")
     parser.add_argument("-o","--outputfile", help="Specify output DX file (default: out.dx) f2-f1",
             required=False, type=str, default="out.dx")
 
     args = parser.parse_args()
     
     if not os.path.isfile(args.inputfile1):
         print("DX File ", args.inputfile1, " does not exist")
         exit(1)

     if not os.path.isfile(args.inputfile2):
         print("DX File ", args.inputfile2, " does not exist")
         exit(1)

     dx1 = Grid(args.inputfile1)

     dx2 = Grid(args.inputfile2)

     g = dx2 - dx1

     g.export(args.outputfile,  type="double")
