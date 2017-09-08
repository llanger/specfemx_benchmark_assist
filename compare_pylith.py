import sys
from gll_library import *
from shape_library import *
from integrate import *
from read_data import *
import settings as g

"""script to compare two specfem3d viscoelastic results, or benchmarking results.
"""

data_file='pylith/pylith_long_fault_diff.csv'
data_len=2884661
connectivity='okada/okada_large_mesh_connectivity'
coordx='okada/okada_large_mesh_coord_x'
coordy='okada/okada_large_mesh_coord_y' 
coordz='okada/okada_large_mesh_coord_z'

#columns containing data in csv file
g.pt_cols=[16,17,18]
g.disp_bench_cols=[6,7,8]
g.disp_sem_cols=[3,4,5]
g.ngllx=3 #change if other number is used
g.ngll=g.ngllx**3

# element size. for use in finding gll pts close to the fault
g.mesh_spacing=

dat,coords,elmt= import_data(data_file,data_len,connectivity,coordx,coordy,coordz)
print "read in all data"
gll_weights,gll_points=gll_quadrature()
print "computed gll quadrature"
dshape_hex8 = dshape_function_hex8(gll_points)
print "computed shape function"
normx,normy,normz,norm=integrate(coords,dat,fault_pts,dshape_hex8,gll_weights,elmt)

f=open('pylith/pylith_long_fault_norms.txt', 'w')
f.write('norm=%3.3f, normx=%3.3f, normy=%3.3f, normz=%3.3f' %(norm,normx,normy,normz))

