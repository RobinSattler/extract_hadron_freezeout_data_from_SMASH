import math
import numpy as np
import sys
import os
import gzip

#time resolution
dt=0.5

#time max
tmax=120

#rapidity limit in absolute value
raplim=1e9

#x,z resolution
dx=0.5
dz=0.5

#histograms go from -xside to xside
xside=20
zside=20

#number of timesteps (automatically set from dt and tmax)
nt=int(math.floor(tmax/dt))

nx=int(math.floor(2*xside/dx))
nz=int(math.floor(2*zside/dz))

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files!=2):
  print ('Syntax: ./compute_dN.py <inputfile> <output_prefix>')
  sys.exit(1)

inputfile=sys.argv[1]
outputprefix=sys.argv[2]

#first, we check that files exist
if(not(os.path.isfile(inputfile))):
  print(inputfile+" does not exist. I quit.\n")
  sys.exit(1)


chemf=0
kinf=1

dN=np.zeros((nt,2),dtype=np.int64)
dN_dxdz=np.zeros((nt,nx,nz,2),dtype=np.int64)

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  pi=gzip.open(inputfile,"r")
else:
  print("Opening file "+inputfile)
  pi=open(inputfile,"r")

hadron_pdg_id=pi.readline().split()[4]
events=int(pi.readline().split()[2])
pi.readline()
pi.readline()
pi.readline()

lines=0
hadrons=np.zeros(2,dtype=np.int64)

for line in pi:
    stuff=line.split() 
    lines=lines+1
    tc,xc,yc,zc,Ec,pxc,pyc,pzc,tk,xk,yk,zk,Ek,pxk,pyk,pzk=np.float64(stuff[:]) 
    try:
        rap_c=0.5*np.log((Ec+pzc)/(Ec-pzc))
    except:
        rap_c=None
    try:
        rap_k=0.5*np.log((Ek+pzk)/(Ek-pzk))
    except:
        rap_k=None
    raps=[rap_c,rap_k]
    t=[tc,tk]
    x=[xc,xk]
    y=[yc,yk]
    z=[zc,zk]
    for q, rap in enumerate(raps):
        if(t[q] < tmax): # the intervals are t_low <= t < t_high
          if(np.abs(rap)<raplim):
           hadrons[q]+=1
           h=int(math.floor(t[q]/dt))
           dN[h,q]+=1
           i=int(math.floor((x[q]+xside)/dx))
           k=int(math.floor((z[q]+zside)/dz))
           if((i>=0) and (i<nx) and (k>=0) and (k<nz)):
                  dN_dxdz[h,i,k,q]=dN_dxdz[h,i,k,q]+1

pi.close()

dN_Ndt=dN/(hadrons*dt) # integers divided by floats are treated like floats
dN_Ndxdz=dN_dxdz/(dx*dz)

print("Total "+hadron_pdg_id+": "+str(lines))
print("Accepted "+hadron_pdg_id+" within rapidity |y|<"+str(raplim)+" at chemical freezeout: "+str(hadrons[chemf]))
print("Accepted "+hadron_pdg_id+" within rapidity |y|<"+str(raplim)+" at kinetic freezeout: "+str(hadrons[kinf]))

#now we print the results into a file
sp="          "
fofl='{:9.5e}'
#foin='{:12d}'
foco='{:5.2f}'
fout=open(outputprefix+"_"+hadron_pdg_id+"_time_distr.txt","w")
fout.write("# rapidity |y|<"+str(raplim)+"\n")
fout.write("# Number of events: "+str(events)+" (averages <> done with respect to the number of events)\n")
fout.write("# Accepted "+hadron_pdg_id+" at chemical freezeout: "+str(hadrons[chemf])+"\n")
fout.write("# Accepted "+hadron_pdg_id+" at kinetic freezeout: "+str(hadrons[kinf])+"\n")
fout.write("#t      dN/(Ndt)(chem f.o.)      <dN/dt>(chem f.o.)     dN/(Ndt)(kin f.o.)      <dN/dt>(kin f.o.)>\n")
for h in range(nt):
    fout.write(foco.format((h+0.5)*dt)+sp+fofl.format(dN_Ndt[h,chemf])+sp+fofl.format(dN[h,chemf]/(dt*events))+sp+fofl.format(dN_Ndt[h,kinf])+sp+fofl.format(dN[h,kinf]/(dt*events))+"\n")
fout.close()

for h in range(nt):
    timestring='{:05.2f}'.format((h+0.5)*dt)
    fout=open(outputprefix+"_"+hadron_pdg_id+"_xz_distr_time_"+timestring+".txt","w")
    fout.write("# rapidity |y|<"+str(raplim)+"\n")
    fout.write("# Number of events: "+str(events)+" (the results are averages across events)\n")
    fout.write("# time: "+timestring+" fm\n")
    fout.write("#x      z      <dN_dxdz>(chem)    <dN_dxdz>(kin)\n")
    for k in range(nz):
        for i in range(nx):
            fout.write(foco.format((i+0.5)*dx)+sp+foco.format((k+0.5)*dz)+sp+fofl.format(dN_dxdz[h,i,k,chemf]/events)+sp+fofl.format(dN_dxdz[h,i,k,kinf]/events)+"\n")
        fout.write("\n")
    fout.close()
