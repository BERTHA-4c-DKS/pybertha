import os
import sys
import argparse

import basicutils

if __name__ == "__main__":

     parser = argparse.ArgumentParser()
     parser.add_argument("-f","--inputfile", help="Specify cube input file", required=True, 
             type=str, default="")
     parser.add_argument("-o","--outputfile", help="Specify DX output file (default: out.dx)",
             required=False, type=str, default="out.dx")
     
     args = parser.parse_args()
     
     if not os.path.isfile(args.inputfile):
         print("CUBE File ", args.inputfile, " does not exist")
         exit(1)

     basicutils.cubetodx (args.inputfile, args.outputfile)


