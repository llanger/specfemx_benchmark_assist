import sys
from gll_library import *
from shape_library import *
from integrate import *
from read_data import *
import settings as g

"""script to compare two specfem3d viscoelastic results, or benchmarking results.
"""

data_file='pylith/pylith_viscoelastic_300yr.csv'
data_len=373979
connectivity='pylith/pylith_specfem_3d_connectivity'
coordx='pylith/pylith_specfem_3d_coord_x'
coordy='pylith/pylith_specfem_3d_coord_y'
coordz='pylith/pylith_specfem_3d_coord_z'
faultfile='pylith/pylith_visco_fault_sources.csv'

#columns containing data in csv file
g.pt_cols=[16,17,18]
g.disp_bench_cols=[10,11,12]
g.disp_sem_cols=[3,4,5]
g.ngllx=3 #change if other number is used
g.ngll=g.ngllx**3

#elmt spacing
g.mesh_spacing=250.0
#number of fault sources
g.npatches=200

dat,coords,elmt,fault_pts= import_data(data_file,data_len,connectivity,coordx,coordy,coordz,faultfile)
print "read in all data"
print 'fault pts', fault_pts
gll_weights,gll_points=gll_quadrature()
print "computed gll quadrature"
dshape_hex8 = dshape_function_hex8(gll_points)
print "computed shape function"
normx,normy,normz,norm=integrate(coords,dat,fault_pts,dshape_hex8,gll_weights,elmt)

f=open('pylith_visco_errors.txt', 'w')
f.write('norm=%3.3f, normx=%3.3f, normy=%3.3f, normz=%3.3f' %(norm,normx,normy,normz))

