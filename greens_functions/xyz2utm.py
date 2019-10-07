import pyproj as pp
import csv 
import sys
""" script to  convert pts from lon/lat to utm coordinates"""

if len(sys.argv) == 4:
    infile_name=str(sys.argv[1])
    outfile_name=str(sys.argv[2])
    xone = float(sys.argv[3])
else:
    raise ValueError("Incorrect number of arguments. Usage: xyz2utm.py <infile_name> <outfile_name> <UTM zone>")    
    exit


putm = pp.Proj(proj='utm', zone=zone, ellps='WGS84')


with open(infile_name, 'rb') as fin, open(outfile_name, 'wb') as fout:
    freader = csv.reader(fin, delimiter='\t')
    fwriter = csv.writer(fout, delimiter = '\t')
    
    for line in freader:
        #convert to UTM 
        line[0], line[1] = putm(float(line[0]),float(line[1]))
        line[2]=float(line[2])
        fwriter.writerow(line)    
