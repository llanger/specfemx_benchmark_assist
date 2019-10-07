import numpy as np
import csv
import settings as g

def import_data(data_file,data_len,connectivity,coordx,coordy,coordz,faultfile):
    with open(data_file,'rb') as fin:
        freader = csv.reader(fin, delimiter=',')
        next(freader, None) #get rid of header line
        i=0
        dat=np.zeros((data_len,9)) #dimension assumed to be 3 for coords + two 3D data cols
        for line in freader:
            for j in range(3): 
                dat[i,j]=line[g.pt_cols[j]]  
                dat[i,j+3]=line[g.disp_bench_cols[j]]
                dat[i,j+6]=line[g.disp_sem_cols[j]]
            i+=1    
    print "finished reading csv data file"        
    with open(faultfile,'rb') as fin:
        freader = csv.reader(fin, delimiter=',')
        next(freader, None) #get rid of header line
        i=0
        fault_pts=np.zeros((g.npatches,3))
        for line in freader:
            fault_pts[i]=line
            i+=1
    print "finished reading fault sources file"    
    with open(connectivity,'rb') as fin:
        g.nelmt=int(fin.readline().strip())
        freader=csv.reader(fin, delimiter=' ')
        elmt=np.zeros((g.nelmt,8),dtype=int)
        i=0
        for line in freader:
            int_line = [int(j) for j in line[0:8]]
            elmt[i]=int_line
            i+=1
    with open(coordx, 'rb') as fin:
        ncoord=int(fin.readline().strip())
        freader=csv.reader(fin,delimiter=' ')
        coords=np.zeros((ncoord,3),dtype=float)
        i=0
        for line in freader:
            coords[i,0]=float(line[0])
            i+=1
    with open(coordy, 'rb') as fin:
        freader=csv.reader(fin,delimiter=' ')
        next(freader,None)
        i=0
        for line in freader:
            coords[i,1]=float(line[0])
            i+=1
    with open(coordz, 'rb') as fin:
        freader=csv.reader(fin,delimiter=' ')
        next(freader,None)
        i=0
        for line in freader:
            coords[i,2]=float(line[0])
            i+=1
    print "finished reading connectivity and coord files"

    return(dat,coords,elmt,fault_pts)    
