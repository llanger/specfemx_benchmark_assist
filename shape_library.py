import numpy as np
import settings as g

def dshape_function_hex8(gll_points):
    """calculate shape function for 8-node element. adapted from shape_library.f90. assumes that ngllx=nglly=ngllz. Does not check that # of nodes is correct (assumes that this has already been checked by specfem)"""

    ngnod=8
    # derivatives of the 3d shape functions
    dshape_hex8 = np.zeros((3,ngnod,g.ngll))

    one_eighth = 0.125
    zerotol=1.0e-12
    # compute the derivatives of 3d shape functions
    igll=0
    for k in range(g.ngllx):
        zetap = 1.0 + gll_points[2,k]
        zetam = 1.0 - gll_points[2,k]
        for j in range(g.ngllx):
            etap = 1.0 + gll_points[1,j]
            etam = 1.0 - gll_points[1,j]
            for i in range(g.ngllx):

                xip = 1.0 + gll_points[0,i]
                xim = 1.0 - gll_points[0,i]

                dshape_hex8[0,0,igll] = - one_eighth*etam*zetam
                dshape_hex8[0,1,igll] = one_eighth*etam*zetam
                dshape_hex8[0,2,igll] = one_eighth*etap*zetam
                dshape_hex8[0,3,igll] = - one_eighth*etap*zetam
                dshape_hex8[0,4,igll] = - one_eighth*etam*zetap
                dshape_hex8[0,5,igll] = one_eighth*etam*zetap
                dshape_hex8[0,6,igll] = one_eighth*etap*zetap
                dshape_hex8[0,7,igll] = - one_eighth*etap*zetap

                dshape_hex8[1,0,igll] = - one_eighth*xim*zetam
                dshape_hex8[1,1,igll] = - one_eighth*xip*zetam
                dshape_hex8[1,2,igll] = one_eighth*xip*zetam
                dshape_hex8[1,3,igll] = one_eighth*xim*zetam
                dshape_hex8[1,4,igll] = - one_eighth*xim*zetap
                dshape_hex8[1,5,igll] = - one_eighth*xip*zetap
                dshape_hex8[1,6,igll] = one_eighth*xip*zetap
                dshape_hex8[1,7,igll] = one_eighth*xim*zetap

                dshape_hex8[2,0,igll] = - one_eighth*xim*etam
                dshape_hex8[2,1,igll] = - one_eighth*xip*etam
                dshape_hex8[2,2,igll] = - one_eighth*xip*etap
                dshape_hex8[2,3,igll] = - one_eighth*xim*etap
                dshape_hex8[2,4,igll] = one_eighth*xim*etam
                dshape_hex8[2,5,igll] = one_eighth*xip*etam
                dshape_hex8[2,6,igll] = one_eighth*xip*etap
                dshape_hex8[2,7,igll] = one_eighth*xim*etap

                igll+=1

    # check the shape functions and their derivatives

    for i in range(g.ngll):
        sum_dshapexi = 0.0
        sum_dshapeeta = 0.0
        sum_dshapezeta = 0.0

        for i_gnod in range(ngnod):
            sum_dshapexi = sum_dshapexi + dshape_hex8[0,i_gnod,i]
            sum_dshapeeta = sum_dshapeeta + dshape_hex8[1,i_gnod,i]
            sum_dshapezeta = sum_dshapezeta + dshape_hex8[2,i_gnod,i]

        # sum of derivative of shape functions should be zero
        if(abs(sum_dshapexi) >  zerotol):
            print'ERROR: derivative xi shape functions!'
            exit 
        if(abs(sum_dshapeeta) >  zerotol):
            print 'ERROR: derivative eta shape functions!'
            exit           
        if(abs(sum_dshapezeta) >  zerotol):
           print 'ERROR: derivative gamma shape functions!'
           exit
    
    return dshape_hex8
