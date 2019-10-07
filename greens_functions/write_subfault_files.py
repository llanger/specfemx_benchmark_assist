import gen_gf_patches
import os
import numpy as np

"""This file is a stripped-down version of the script that I used to write the moment tensor solution files required for SPECFEMX. Unfortunately I can't share the part of the code that calculates the geometry for each patch. However, if you have your own code to produce the geometry this will write the output in the proper format.""" 

npatch=20 #number of patches per dimension in the subfault. This will generate 20 * 20 = 400 patches. Use more or less for bigger / smaller subfaults.
subfaults = np.genfromtxt('yourgeometryfile.csv', delimiter = ',')

# origin of x,y coordinate system for your mesh
xo=0.0
yo=0.0

#topographic elevation of the lowest point on the hinge line. The fault should begin here if it intersects the surface
zo=244.0

i=0
for p in ptc.patch:
    xc, yc, zc, width, length, vs, vd, vn = subfaults[i] 
    for j in range(2):
        if j==0:
            f=open("patches/patch%s_ds" % i, "w")
        elif j==1:
            f=open("patches/patch%s_ss" % i, "w")
        f.write("# slip information file for patch %i \n" % i)
        f.write("fault_center = %3.3f %3.3f %3.3f \n" % (xc-xo,yc-yo,zc+zo))
        f.write("fault_strike = ")
        np.savetxt(f,vs[None],fmt='%3.3f')
        f.write("fault_normal = ")
        np.savetxt(f,vn[None],fmt='%3.3f')
        f.write("fault_slip = ")
        if j==0:
            np.savetxt(f,-1*vd[None],fmt='%3.3f')
        elif j ==1:
            np.savetxt(f,vs[None],fmt='%3.3f')
        f.write("fault_length = %3.3f \n" % length)
        f.write("fault_width = %3.3f \n" % width)
        f.write("npatch_length = %3.3f \n" % npatch)
        f.write("npatch_width = %3.3f \n" % npatch)
        f.close()
    i+=1
