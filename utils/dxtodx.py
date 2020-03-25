import os
import sys
import argparse

sys.path.insert(0, "/home/redo/Sources/GridDataFormats/")

from gridData import Grid

import numpy as np

import basicutils

if __name__ == "__main__":

     parser = argparse.ArgumentParser()
     parser.add_argument("-f","--inputfile", help="Specify DX input file", required=True, 
             type=str, default="")
     parser.add_argument("-o","--outputfile", help="Specify output DX file (default: out.dx)",
             required=False, type=str, default="out.dx")
 
     args = parser.parse_args()
     
     if not os.path.isfile(args.inputfile):
         print("DX File ", args.inputfile, " does not exist")
         exit(1)

     fp = open(args.inputfile)

     headerdone = False
     header = []
     values = []
     tail = []
     for line in fp:
         if line != "\n":
             if line.find("object") >= 0 or line.find("orig") >= 0 or \
                     line.find("delta") >= 0 or line.find("component") >= 0 \
                     or line.find("end") >= 0 or line.find("#") >= 0 or \
                     line.find("attribute") >= 0:
                 if headerdone:
                     tail.append(line.lstrip().rstrip())
                 else:
                     header.append(line.lstrip().rstrip())
             else:
                 headerdone = True
                 values.append(np.float64(line))
     
     fp.close()

     fp = open(args.outputfile, "w")

     fp.write("# Original file " + args.inputfile + "\n")
     fp.write("# converted\n")
     for l in header:
         if l != "end":
            fp.write(l+"\n")
     counter = 0
     for e in values:
         fp.write("%20.15f "%(e))
         counter += 1
         if counter == 3:
             fp.write("\n")
             counter = 0

     if counter < 3:
         fp.write("\n")

     for l in tail:
         if l != "end":
            fp.write(l+"\n")

     fp.close()

     g = Grid(args.outputfile)

     g.export(args.outputfile,  type="double")
