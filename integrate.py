import numpy as np
import numpy_indexed as npi
import settings as g
from scipy.spatial import distance
from math import sqrt

def integrate(coords,data,fault_pts,dshape_hex8,gll_weights,elmt):
    """function to perform the integration for metric computation"""
    norm=0.0
    normx=0.0
    normy=0.0
    normz=0.0
    div=0.0 #normalizing factor to divide by
    divx=0.
    divy=0.
    divz=0.

    eps=1.0*g.mesh_spacing/(g.ngllx-1.)
    print 'eps=', eps
    f=open('eliminated_coords.vtk','w')

    #create integer versions of arrays to use in pulling out gll pts for each element
    data_round=np.rint(data)
    dati=data_round.astype(int)
    coord_round=np.rint(coords)
    coordi=coord_round.astype(int)

    #remove duplicates from data array
    dat_struc=np.ascontiguousarray(dati).view(np.dtype((np.void,dati.dtype.itemsize *dati.shape[1])))
    _,idx=np.unique(dat_struc,return_index=True)
    datu=dati[idx]
    data_unique=data[idx]

    for i_elmt in range(g.nelmt):
        #pull out geometric coordinates for this element
        elmt_coord_id=[j-1 for j in elmt[i_elmt]]
        elmt_coord=coordi[elmt_coord_id]

        #find corresponding gll pts for this element
        xmin=min(elmt_coord[:,0]);xmax=max(elmt_coord[:,0])
        ymin=min(elmt_coord[:,1]);ymax=max(elmt_coord[:,1])
        zmin=min(elmt_coord[:,2]);zmax=max(elmt_coord[:,2])
        gll_coord_id=np.nonzero((datu[:,0]>=xmin) & (datu[:,0]<=xmax) & (datu[:,1]>=ymin) & (datu[:,1]<=ymax) & (datu[:,2]>=zmin) & (datu[:,2]<=zmax))
        elmt_data=data_unique[gll_coord_id]
        if len(gll_coord_id[0]) != g.ngll:
            print "elmt=", elmt_coord_id
            print xmin,xmax,ymin,ymax,zmin,zmax
            print 'elmt_data=', elmt_data
            print "gll pts found=", len(gll_coord_id[0])
            raise ValueError("incorrect number of gll points found in element!")
            exit

        #sort the gll coords so they correspond the order of the arrays giving the weights and shape function
        dat_sorted=elmt_data[npi.argsort((elmt_data[:,0], elmt_data[:,1],elmt_data[:,2]))]
        func=dat_sorted[:,3:]

        #if any gll pt is too close to fault, remove the element from the integration
        dist=distance.cdist(fault_pts,dat_sorted[:,0:3],'euclidean')
        if (dist<eps).any():
            print "eliminated element #", i_elmt
            np.savetxt(f,dat_sorted[:,0:3],fmt='%3.3f')
            continue

        for i_gll in range(g.ngll):

            #compute jacobian, its derivative and inverse
            jac=np.matmul(dshape_hex8[:,:,i_gll],elmt_coord)
            det_jac=np.linalg.det(jac)

            #perform the integration
            norm=norm+det_jac*gll_weights[i_gll]*np.dot((func[i_gll,3:6]-func[i_gll,0:3]),(func[i_gll,3:6]-func[i_gll,0:3]))
            div=div+det_jac*gll_weights[i_gll]*np.dot(func[i_gll,3:6],func[i_gll,3:6])
            normx=normx+det_jac*gll_weights[i_gll]*(func[i_gll,3]-func[i_gll,0])**2
            divx=divx+det_jac*gll_weights[i_gll]*(func[i_gll,3])**2
            normy=normy+det_jac*gll_weights[i_gll]*(func[i_gll,4]-func[i_gll,1])**2
            divy=divy+det_jac*gll_weights[i_gll]*(func[i_gll,4])**2
            normz=normz+det_jac*gll_weights[i_gll]*(func[i_gll,5]-func[i_gll,2])**2
            divz=divz+det_jac*gll_weights[i_gll]*(func[i_gll,5])**2

    norm_finalx=sqrt(normx/divx)
    norm_finaly=sqrt(normy/divy)
    norm_finalz=sqrt(normz/divz)
    norm_final=sqrt(norm/div)

    f.close()

    return norm_finalx, norm_finaly, norm_finalz,norm_final
            


