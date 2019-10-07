import numpy as np
import pyproj as pp

"""read in GPS and InSAR data files and write a file containing the UTM coordinates of the data points in .vtk format"""

# Load data from CSV
#GPS data
dat_gps = np.genfromtxt('your_gps_data_file.txt',skip_header=0)
gps_loc=dat_gps[:,1:3] #change this if your file format is different

#INSAR  data
dat_insar = np.genfromtxt('your_insar_data.txt',skip_header=1)
insar_loc=dat_insar[:,3:5]

#remove coords outside domain
def remove_outside_coords(loc):
    nulls=[]
    for i in range(len(loc)):
        if (loc[i,0] < 83.5 or loc[i,0] > 87.47 or loc[i,1] < 26.6 or loc[i,1] > 29.2): #replace these numbers with lat/lon limits of your mesh domain
            nulls.append(i)
    loc=np.delete(loc,nulls,axis=0)
    print nulls
    return loc

gps_loc=remove_outside_coords(gps_loc)
insar_loc=remove_outside_coords(insar_loc)

print "number of data points for GPS/ insar:", len(gps_loc),len(insar_loc)

#combine all of the lat/lon coordinates into a single array
locs=np.vstack((gps_loc,insar_loc))
#initialize empty array to hold UTM coordinates
coords=np.zeros((len(locs),3))

# origin of x,y coordinate system for the mesh. should only be nonzero if you've shifted relative to standard UTM
xo=0.0
yo=0.0

#convert to utm
putm = pp.Proj(proj='utm', zone=45, ellps='WGS84')#change to proper UTM zone
for i in range(len(locs)):
    x,y=putm(locs[i,0],locs[i,1])
    coords[i,0]=x-xo
    coords[i,1]=y-yo

pts=len(coords)
#generate array for vtk file
num=np.ones((pts,2),dtype=int)
for i in range(pts):
        num[i,1]=i+1

#write to vtk file
f=open('station_coords.vtk', 'w')
f.write("# vtk DataFile Version 2.0\n")
f.write("stations\n")
f.write("ASCII\n")
f.write("DATASET POLYDATA\n")
f.write("POINTS %d float\n" %(pts))
np.savetxt(f,coords,fmt='%3.3f')
f.write("\n")
f.write("VERTICES %d %d\n" %(pts,pts*2))
np.savetxt(f,num.astype(int),fmt='%i')
f.close()

