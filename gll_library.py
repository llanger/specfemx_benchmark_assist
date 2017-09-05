from scipy.special import gamma
from math import atan,cos
import numpy
import settings as g
""" functions to calculate gll pts, weights, lagrange polynomials, etc. Modified from gll_library.f90 in specfem3d_viscoelastic, written by HN Gharti"""

def jacobf(n,alp,bet,x):

    apb  = alp+bet
    poly = 1.0
    pder = 0.0
    psave = 0.0
    pdsave = 0.0

    if (n == 0): return

    polyl = poly
    pderl = pder
    poly  = (alp-bet+(apb+2.0)*x)/2.0
    pder  = (apb+2.0)/2.0
    if (n == 1): return

    for k in range(2,n+1):
        dk = k
        a1 = 2.0*dk*(dk+apb)*(2.0*dk+apb-2.0)
        a2 = (2.0*dk+apb-1.0)*(alp**2-bet**2)
        b3 = (2.0*dk+apb-2.0)
        a3 = b3*(b3+1.0)*(b3+2.0)
        a4 = 2.0*(dk+alp-1.0)*(dk+bet-1.0)*(2.0*dk+apb)
        polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
        pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
        psave  = polyl
        pdsave = pderl
        polyl  = poly
        poly   = polyn
        pderl  = pder
        pder   = pdern

    polym1 = polyl
    pderm1 = pderl
    polym2 = psave
    pderm2 = pdsave

    return poly,pder,pderm1

def endw1(n,alpha,beta):

    f3 = 0.0
    apb   = alpha+beta
    if (n == 0):
        return 0.0

    f1   = gamma(alpha+2.0)*gamma(beta+1.0)/gamma(apb+3.0)
    f1   = f1*(apb+2.0)*2.0**(apb+2.0)/2.0
    if (n == 1):
        return f1
   
    fint1 = gamma(alpha+2.0)*gamma(beta+1.0)/gamma(apb+3.0)
    fint1 = fint1*2.0**(apb+2.0)
    fint2 = gamma(alpha+2.0)*gamma(beta+2.0)/gamma(apb+4.0)
    fint2 = fint2*2.0**(apb+3.0)
    f2    = (-2.0*(beta+2.0)*fint1 + (apb+4.0)*fint2) * (apb+3.0)/4.0
    if (n == 2):
        return f2
   
    for i in range(3,n+1):
        di   = i-1
        abn  = alpha+beta+di
        abnn = abn+di
        a1   = -(2.0*(di+alpha)*(di+beta))/(abn*abnn*(abnn+1.0))
        a2   =  (2.0*(alpha-beta))/(abnn*(abnn+2.0))
        a3   =  (2.0*(abn+1.0))/((abnn+2.0)*(abnn+1.0))
        f3   =  -(a2*f2+a1*f1)/a3
        f1   = f2
        f2   = f3
  
    return f3

def endw2(n,alpha,beta):

    apb   = alpha+beta
    f3 = 0.0
    if (n == 0):
        return 0.0
    
    f1   = gamma(alpha+1.0)*gamma(beta+2.0)/gamma(apb+3.0)
    f1   = f1*(apb+2.0)*2.0**(apb+2.0)/2.0
    if (n == 1):
        return f1
    
    fint1 = gamma(alpha+1.0)*gamma(beta+2.0)/gamma(apb+3.0)
    fint1 = fint1*2.0**(apb+2.0)
    fint2 = gamma(alpha+2.0)*gamma(beta+2.0)/gamma(apb+4.0)
    fint2 = fint2*2.0**(apb+3.0)
    f2    = (2.0*(alpha+2.0)*fint1 - (apb+4.0)*fint2) * (apb+3.0)/4.0
    if (n == 2):
        return f2
    
    for i in range(3,n+1):
        di   = i-1
        abn  = alpha+beta+di
        abnn = abn+di
        a1   =  -(2.0*(di+alpha)*(di+beta))/(abn*abnn*(abnn+1.0))
        a2   =  (2.0*(alpha-beta))/(abnn*(abnn+2.0))
        a3   =  (2.0*(abn+1.0))/((abnn+2.0)*(abnn+1.0))
        f3   =  -(a2*f2+a1*f1)/a3
        f1   = f2
        f2   = f3

    return f3

def  jacg(np,alpha,beta):
    
    pm1 = 0.0
    pm2 = 0.0
    pdm1 = 0.0
    pdm2 = 0.0
    eps=1.0e-12
    xlast = 0.0
    n   = np-1
    dth = 4.0*atan(1.0)/(2.0*n+2.0)
    p = 0.0
    pd = 0.0
    jmin = 0
    xjac=numpy.zeros(np)
    for j in range(1,np+1):
        if(j == 1):
            x = cos((2.0*(j-1.0)+1.0)*dth)
        else:
            x1 = cos((2.0*(j-1.0)+1.0)*dth)
            x2 = xlast
            x  = (x1+x2)/2.0
        for k in range(20):
            p,pd,pdm1=jacobf(np,alpha,beta,x)
            recsum = 0.0
            jm = j-1
            for i in range(1,jm+1):
                recsum = recsum+1.0/(x-xjac[np-i])
            delx = -p/(pd-recsum*p)
            x    = x+delx
            if(abs(delx) < eps): break
        xjac[np-j] = x
        xlast        = x
    for i in range(1,np+1):
        xmin = 2.0
        for j in range(i,np+1):
            if(xjac[j-1] < xmin):
                xmin = xjac[j-1]
                jmin = j 
        if(jmin != i):
            swap = xjac[i-1]
            xjac[i-1] = xjac[jmin-1]
            xjac[jmin-1] = swap

    return xjac

def lagrange1d(nenode,xi):

    #initialize arrays
    phi=numpy.zeros(nenode);dphi_dxi=numpy.zeros(nenode)
    xii=numpy.zeros(nenode);term=numpy.zeros(nenode);
    dterm=numpy.zeros(nenode);sum_term=numpy.ones(nenode)

    # compute natural coordnates
    dx=2.0/(nenode-1.)# length = 2.0 as xi is taken -1 to +1
    for i in range(nenode):
        # coordinates when origin is in the left
        xii[i]=(i-1.)*dx

    # origin is tranformed to mid point
    xii=xii-1.0

    for i in range(nenode):
        k=0
        phi[i]=1.0
        for j in range(nenode):
            if(j!=i):
                term[k]=(xi-xii[j])/(xii[i]-xii[j])
                dterm[k]=1.0/(xii[i]-xii[j]) # derivative of the term wrt xi
                phi[i]=phi[i]*(xi-xii[j])/(xii[i]-xii[j])
                k=k+1
        for j in range(nenode-1):
            for k in range(nenode-1):
                if(k==j):
                    sum_term[j]=sum_term[j]*dterm[k]
                else:
                    sum_term[j]=sum_term[j]*term[k]
        dphi_dxi[i]=0.0
        for j in range(nenode-1):
            dphi_dxi[i]=dphi_dxi[i]+sum_term[j]

    return phi,dphi_dxi

def pnormj(n,alpha,beta):
    const = alpha+beta+1.0

    if (n <= 1):
        prod   = gamma(n+alpha)*gamma(n+beta)
        prod   = prod/(gamma(n)*gamma(n+alpha+beta))
        pnormj = prod * 2.0**const/(2.0*n+const)
        return pnormj

    prod  = gamma(alpha+1.0)*gamma(beta+1.0)
    prod  = prod/(2.0*(1.0+const)*gamma(const+1.0))
    prod  = prod*(1.0+alpha)*(2.0+alpha)
    prod  = prod*(1.0+beta)*(2.0+beta)

    for i in range(3,n+1):
        frac  = (i+alpha)*(i+beta)/(i*(i+alpha+beta))
        prod  = prod*frac

    pnormj = prod * 2.0**const/(2.0*n+const)
    return pnormj


def gll_quadrature():
    """ generate gll points and weights for the given gll coordinates. assumes all ngllx, nglly, ngllz are all the same -- so only ngllx is specified. also assumes that ngllx is at least 3. Adapted from gll_library.f90"""

    alpha=0.0
    beta=0.0
    gll_weights = numpy.zeros(g.ngll)
    gll_points=numpy.zeros((3,g.ngll))
    gllp=numpy.zeros(g.ngllx)
    gllw=numpy.zeros(g.ngllx)
    
    alpg=alpha+1.0; betg=beta+1.0
    n=g.ngllx-1
    nm1=g.ngllx-2
    apb=alpg+betg
    
    if(nm1==1): 
        gllp[1]=(betg-alpg)/(apb+2.0)
        gllw[1]=gamma(alpg+1.0)*gamma(betg+1.0)/gamma(apb+2.0) *2.0**(apb+1.0)
    else: 
        gllp[1:n]=jacg(nm1,alpg,betg)
        np1   = nm1 
        np2   = nm1+1
        fac1  = np1+alpg+betg+1.0
        fac2  = fac1+np1
        fac3  = fac2+1.0
        fnorm = pnormj(np1,alpg,betg)
        rcoef = (fnorm*fac2*fac3)/(2.0*fac1*np2)
        for i in range(1,n):
            p,pd,pdm1=jacobf(np2,alpg,betg,gllp[i])
            gllw[i] = -rcoef/(p*pdm1)

    gllp[0]=-1.0
    gllp[g.ngllx-1]=1.0

    for  i in range(1,n):
      gllw[i] = gllw[i]/(1.0-gllp[i]**2)

    p,pd,pdm1=jacobf(n,alpha,beta,gllp[0])
    gllw[0]  = endw1(n,alpha,beta)/(2.0*pd)
    p,pd,pdm1=jacobf(n,alpha,beta,gllp[g.ngllx-1])
    gllw[g.ngllx-1] = endw2(n,alpha,beta)/(2.0*pd)

    l=0
    for k in range(g.ngllx):
        for j in range(g.ngllx):
            for i in range(g.ngllx):
                # integration points -- assumes gllpx, gllpy, gllpz are identical
                gll_points[0,l]=gllp[i]
                gll_points[1,l]=gllp[j]
                gll_points[2,l]=gllp[k]

                # integration weights
                gll_weights[l]=gllw[i]*gllw[j]*gllw[k]
                l+=1

    #compute 1d lagrange polynomials & derivatives            
    """    dlagrange_gll=numpy.zeros((3,g.ngll,g.ngll))
    for i_gll in range(g.ngll):

        lagrange_x,lagrange_dx=lagrange1d(g.ngllx,gll_points[0,i_gll])
        lagrange_y,lagrange_dy=lagrange1d(g.ngllx,gll_points[1,i_gll])
        lagrange_z,lagrange_dz=lagrange1d(g.ngllx,gll_points[2,i_gll])

        n=0 
        for k in range(g.ngllx):
            for j in range(g.ngllx):
                for i in range(g.ngllx):
                   #lagrange_gll(ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_z(k)
                   dlagrange_gll[0,i_gll,n]=lagrange_dx[i]*lagrange_y[j]*lagrange_z[k]
                   dlagrange_gll[1,i_gll,n]=lagrange_x[i]*lagrange_dy[j]*lagrange_z[k]
                   dlagrange_gll[2,i_gll,n]=lagrange_x[i]*lagrange_y[j]*lagrange_dz[k]
                   n=n+1
    """
    return gll_weights,gll_points
    
