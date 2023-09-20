#!/usr/bin/python3

# get_fo_data.py 05/04/2022
# 
# it reads the hadron data in the smash collision files and extracts
# the position and the momenta of selected hadrons at chemical and kinetic freezeout

import math
import numpy as np
import os
import sys

# tuple with the PDG IDs of the hadron to search
hadron_tuple = (3122,3212)

# we ignore time differences smaller than this one
tiny_dt=1.e-12

# format for quantities in output file
ff='{:14.10e}'
sp="    "

# if we want to print debugging messages or not (0=none,1=medium,2=all)
#verbose = 2

# indexes when reading entry lines of data collision file
i_t = 0
i_x = 1
i_y = 2
i_z = 3
i_zR = i_z + 1 # for ending ranges

i_E = 5
i_px = 6
i_py = 7
i_pz = 8
i_pzR = i_pz + 1 # for ending ranges

i_PDG_id = 9
i_simulation_id = 10
i_form_time = 13

particle_data = {}

class fo_data:
    def __init__(self,param):
        self.pdg_id = param
        self.form_pos = np.zeros(4,dtype=np.float64)
        self.form_mom = np.zeros(4,dtype=np.float64)
        self.dec_pos = np.zeros(4,dtype=np.float64)
        self.dec_mom = np.zeros(4,dtype=np.float64)

def process_outgoing(dataline):
    particle_simulation_ID = dataline[i_simulation_id]
    position = np.float64(dataline[i_t:i_zR])
    momentum = np.float64(dataline[i_E:i_pzR])
    form_time = float(dataline[i_form_time])
    dt = form_time - position[0]
    if dt > tiny_dt:
        En = np.float64(dataline[i_E])
        vel = momentum[1:4]/En
        position[1:4] = position[1:4] + vel*dt
        position[0] = form_time

    if ( not particle_simulation_ID in particle_data ):
        #if (verbose > 1):
        #    print("Adding new hadron with id "+particle_simulation_ID)
        particle_data[particle_simulation_ID] = fo_data(int(dataline[i_PDG_id]))
        #particle_data[particle_simulation_ID].pdg_id = int(dataline[i_PDG_id])
        particle_data[particle_simulation_ID].form_mom[:] = momentum[:]
        particle_data[particle_simulation_ID].form_pos[:] = position[:]

    #if (verbose > 1):
    #    print("Updating decoupling info of hadron with id "+particle_simulation_ID)
    particle_data[particle_simulation_ID].dec_pos[:] = position[:]
    particle_data[particle_simulation_ID].dec_mom[:] = momentum[:]
    #if (verbose > 1):
    #    print("Four momentum at decoupling: "+str(particle_data[particle_simulation_ID].dec_mom[:]))


def extract_data(inputfile,hadron_list):
    try:
         infile=open(inputfile,"r")
    except OSError:
         print("Could not open/read file: ", inputfile)
         sys.exit(1)
   
    n_events = 0
    results = []
    for line in infile:
        stuff=line.split()
        if(len(stuff)==0):
            continue
        if stuff[1] == "interaction": 
            in_part = int(stuff[3])
            out_part = int(stuff[5])
            potentially_to_remove=[]
    #        if (verbose > 1):
    #            print("Interaction: in "+stuff[3]+" out "+stuff[5])
            for iitt in range(in_part):
                stuff=infile.readline().split()
    #            if (verbose > 1):
    #                print("In "+stuff[i_PDG_id])
                if int(stuff[i_PDG_id]) in hadron_list:
    #                if (verbose > 1):
    #                    print("Observed "+stuff[i_PDG_id]+" with sim id: "+stuff[i_simulation_id])
                    potentially_to_remove.append(stuff[i_simulation_id])
            for iitt in range(out_part):
                stuff=infile.readline().split()
    #            if (verbose > 1):
    #                print("Out "+stuff[i_PDG_id])
                if int(stuff[i_PDG_id]) in hadron_list:
                    process_outgoing(stuff)
    #                if (verbose > 0):
    #                    print("Added "+stuff[i_PDG_id]+" with sim id: "+stuff[i_simulation_id])
                    if stuff[i_simulation_id] in potentially_to_remove:
    #                    if (verbose > 0):
    #                        print("Keeping "+stuff[i_PDG_id]+" with sim id: "+stuff[i_simulation_id])
                        potentially_to_remove.remove(stuff[i_simulation_id])
            if (len(potentially_to_remove)>0):
                for entry in potentially_to_remove:
    #                if (verbose > 0):
    #                    print("DBG Removing "+str(entry))
                    del particle_data[entry]
    #        if (verbose > 1):
    #            print("DBG: dictionary length: "+str(len(particle_data)))
    #            for key in particle_data.keys():
    #                print("DBG dict: key: "+str(key)+" particle: "+str(particle_data[key].pdg_id)+" momentum: "+str(particle_data[key].form_mom[:]))
        if stuff[1] == "event": 
            n_events = n_events + 1
            for key, value in particle_data.items():
    #            if (verbose > 0):
    #                print("DBG event: "+str(n_events)+" key: "+str(key)+" particle: "+str(value.pdg_id)+" momentum: "+str(value.form_mom[:]))
                results.append(value)
            particle_data.clear()
    return n_events, results

if (__name__ == "__main__" ):
    if (len(sys.argv)<3):
        print ('Syntax: ./get_fo_data.py <output file suffix> <smash collision file 1> [smash collision file 2] ...')
        sys.exit(1)
    else:
        total_events=0
        total_results=[]
        outputfiles=[]
        for h in hadron_tuple:
            outputfiles.append(str(h)+"_"+sys.argv[1])
        for hf in outputfiles:
            if (os.path.exists(hf)):
                print("Output file "+hf+" already exists, I will not overwrite it. I stop here.")
                sys.exit(2)
        
        for afl in sys.argv[2:]:
            new_events, new_results = extract_data(afl, hadron_tuple)
            total_events = total_events + new_events
            total_results.extend(new_results)

        outf={}
        for oindx, hf in enumerate(outputfiles):
            h_tuple_id=hadron_tuple[oindx]
            outf[h_tuple_id]=open(hf,"w")
            outf[h_tuple_id].write("# hadron PDG ID: "+str(hadron_tuple[oindx])+"\n")
            outf[h_tuple_id].write("# events: "+str(total_events)+"\n")
            outf[h_tuple_id].write("# columns 1-4 time and position at formation, columns 5-8 four momenta at formation\n")
            outf[h_tuple_id].write("# columns 9-12 time and position at last collision, columns 13-15 four momenta after last collision\n")
            outf[h_tuple_id].write("# t x y z E px py pz t' x' y' z' E' px' py' pz'\n")
        
        for entry in total_results:
            for ll in range(0,4):
                outf[entry.pdg_id].write(ff.format(entry.form_pos[ll])+sp)
            for ll in range(0,4):
                outf[entry.pdg_id].write(ff.format(entry.form_mom[ll])+sp)
            for ll in range(0,4):
                outf[entry.pdg_id].write(ff.format(entry.dec_pos[ll])+sp)
            for ll in range(0,4):
                outf[entry.pdg_id].write(ff.format(entry.dec_mom[ll])+sp)
            outf[entry.pdg_id].write("\n")


        for hf in outf.values():
            hf.close()
