import sys
from gll_library import *
from shape_library import *
from integrate import *
from read_data import *
import settings as g
import numpy as np

np.set_printoptions(threshold=np.inf)

"""script to compare two specfem3d viscoelastic results, or benchmarking results.
"""

data_file='okada/hn_okada_cmt.csv'
data_len=433593
connectivity='okada/vertical_fault_connectivity'
coordx='okada/vertical_fault_coord_x'
coordy='okada/vertical_fault_coord_y'
coordz='okada/vertical_fault_coord_z'
faultfile='okada/hn_okada_cmt_sources.csv'

#columns containing data in csv file
g.pt_cols=[9,10,11]
g.disp_bench_cols=[6,7,8]
g.disp_sem_cols=[0,1,2]
g.ngllx=3 #change if other number is used
g.ngll=g.ngllx**3

#elmt spacing
g.mesh_spacing=1000.0
#number of fault sources
g.npatches=200

dat,coords,elmt,fault_pts= import_data(data_file,data_len,connectivity,coordx,coordy,coordz,faultfile)
print "read in all data"
gll_weights,gll_points=gll_quadrature()
print "computed gll quadrature"
dshape_hex8 = dshape_function_hex8(gll_points)
print "computed shape function"
normx,normy,normz,norm=integrate(coords,dat,fault_pts,dshape_hex8,gll_weights,elmt)

f=open('okada/hn_okada_cmt_errors.txt', 'w')
f.write('norm=%3.3f, normx=%3.3f, normy=%3.3f, normz=%3.3f' %(norm,normx,normy,normz))

