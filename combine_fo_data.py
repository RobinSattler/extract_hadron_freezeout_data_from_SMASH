#!/usr/bin/python3

# combine_fo_data.py 09/04/2022
# 
# it combines the output files produced by get_fo_data.py

import math
import numpy as np
import os
import sys

if len(sys.argv) < 4:
    print("Syntax: python3 combine_fo_data.py <outputfile> <inputfile1> <inputfile2> ... <inputfileN>")
    sys.exit(1)

if (os.path.exists(sys.argv[1])):
    print("Output file "+sys.argv[1]+" already exists. I will no overwrite it.")
    sys.exit(1)



total_events = 0
for ii, infx in enumerate(sys.argv[2:]):
    if (not os.path.exists(infx)):
        print("Warning, input file "+infx+" does not exist.")
    else:
        infile=open(infx,"r")
        pdg_id = infile.readline().split()[4]
        total_events = total_events + int(infile.readline().split()[2])
        infile.close()
        if ii == 0:
            pdg_id_ref = pdg_id
        else:
            if (pdg_id != pdg_id_ref):
                print("Error, input file "+infx+" contains data about "+pdg_id+" instead of "+pdg_id_ref+"!!! I quit.")
                sys.exit(2)

with open(sys.argv[1],"w") as outf:
    outf.write("# hadron PDG ID: "+pdg_id_ref+"\n")
    outf.write("# events: "+str(total_events)+"\n")
    outf.write("# columns 1-4 time and position at formation, columns 5-8 four momenta at formation\n")
    outf.write("# columns 9-12 time and position at last collision, columns 13-15 four momenta after last collision\n")
    outf.write("# t x y z E px py pz t' x' y' z' E' px' py' pz'\n")

    for inputfile in sys.argv[2:]:
        infile=open(inputfile,"r")
        for i in range(5):
            infile.readline()
        for line in infile:
            outf.write(line)
        infile.close()
