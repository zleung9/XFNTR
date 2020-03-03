import sys
import math as m
import numpy as np
from numba import jit
from numba.typed import List
import matplotlib.pyplot as plt
import periodictable
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from scipy.special import erf

e_radius = periodictable.constants.electron_radius*1e10
energy = 20 # keV


#### data generating ####

def sldCalFun(d,rho,sigma,x):

    N = len(rho) - 1 # number of internal interfaces
    z = [0] + [sum(d[:j+1]) for j in range(len(d))] # height of each interface
    sld_tot = np.zeros(x.shape)
    
    # sld calculation model in Wei's paper 
    for i in range(N):
        erfx = (x-z[i])/sigma[i]/np.sqrt(2)
        sld = 0.5 * (erf(erfx) * (rho[i+1]-rho[i]))
        sld_tot = sld_tot + sld
    sld_tot = sld_tot + 0.5 * (rho[0]+rho[N])

    return sld_tot

def refCalFun(d,rho,mu,sigma,x,rrf=False):
    
    if sigma[0] <= 0: sigma[0] = 1e-5 # eliminate zero
    sigma = sigma[0] * np.ones(len(sigma)) # fixed sigma

    erad = e_radius
    slab=0.25

    x_ = List()
    for xx in x: x_.append(xx)

    k0=2*np.pi*float(energy)/12.3984 # wave vector
    lamda=2*np.pi/k0
    theta=x/2/k0   # convert q to theta
   
    # total length of inner slabs plus 4 times roughness for both sides
    length = np.sum(d) + 4* (sigma[0]+sigma[-1])
    steps=int(length/slab) # each sliced box has thickness of ~ 0.25 A
    xsld=np.linspace(-4*sigma[0],np.sum(d)+4*sigma[-1],steps) # get the z-axis for sld

    sd=length/steps # thickness for each slab
    intrho=sldCalFun(d,rho,sigma,xsld) # rho for all the steps
    intmu=sldCalFun(d,mu,sigma,xsld) # mu for all the steps
    

    # was....
    d, sdel, sbet = List(), List(), List()
    sdel.append(float(rho[0])) # delta for the top phase
    sbet.append(float(mu[0]/2/k0/1e8))        # beta for the top phase
    for rho_ in intrho: # add rho for the interface
        sdel.append(float(rho_))
    sdel.append(float(rho[-1]))  # delta for the bottom phase
    # add delta for the interface
    for mu_ in intmu: # add beta for the interface
        sbet.append(float(mu_/2/k0/1e8))
    sbet.append(float(mu[-1]/2/k0/1e8))    # beta for the bottom phase         
    for i in sdel: d.append(slab)
    ref,refr=parratt(x_, lamda, d, sdel, sbet)

    d_, sdel_, sbet_ = List(), List(), List()
    for i in [0.0,1.0]: d_.append(i)
    for i in [sdel[0],sdel[-1]]: sdel_.append(i)
    for i in [sbet[0],sbet[-1]]: sbet_.append(i)
    frsnll,frsnl1=parratt(x_, lamda,d_, sdel_, sbet_) 

    if rrf == True:
        return ref/frsnll
    else:
        return ref

def ref2min(params,x,y,yerr,fit=True):

    ndata = len(x)

    # allocate parameters to different lists
    sigma_t, rho_b, mu = [], [], []
    d, sigma, rho, qoff = [], [], [], []

    par = params.valuesdict()
    for p in par:
        if p.startswith('sigma'):
            if p.endswith('0'):
                sigma_t.append(par[p]); continue
            else:
                sigma.append(par[p]); continue
        if p.startswith('d'):
            d.append(par[p]); continue
        if p.startswith('mu'):
            mu.append(par[p]); continue
        if p.startswith('qoff'):
            qoff.append(par[p]); continue
        if p.startswith('rho'):
            if '_b' in p:
                rho_b.append(par[p]); continue
            else:
                rho.append(par[p]); continue

    if fit==True: # return residual: 1D array
        residual = np.array([])
        for i in range(ndata):
            yy = refCalFun(d,
                           rho+[rho_b[i]],
                           mu,
                           [sigma[i]]+sigma_t,
                           x[i])
            residual = np.append(residual, (yy-y[i])/yerr[i])             
        return residual
    else: # return model: (y1,y2,y3,y4,y5)
        model = []
        for i in range(ndata):
            yy = refCalFun(d,
                           rho + [rho_b[i]],
                           mu,
                           [sigma[i]] + sigma_t,
                           x[i])
            model.append(yy)
        return tuple(model)
   
def iterCallBack(params,iteration,resid,x,y,err,fit=True):
    if fit== False: return None
    m = sum([params[p].vary for p in params]) # number of parameters
    n = sum([len(xx) for xx in x]) # number of points
    df = n - m # degree of freedom
    redchisq = sum(resid**2) / df
    if (iteration<=10) or (iteration%10==0):  #display reduced chisq every 10 iteration
        print(iteration,redchisq)

def initParameters(name_list,para_list):
    
    if len(name_list)!=len(para_list):
        print(" Does name_list and value_list match, please check ")
        return
    params = Parameters()
    for i,name in enumerate(name_list):
        p = para_list[i]
        params.add(name,value=p[0],vary=p[1],min=p[2],max=p[3])
    return params   

def updateParameters(mrefpar,refparaname,refpara):
    '''Return the parameter name for the items selected in the table.'''
    ui = mrefpar  # it only accepts ref_multiFit_par.ui as ui
    
    parameter_list = refpara
    name_list = refparaname
    index_dict = name2index(name_list,reverse=True)

    # update the value of each cell in the table
    for index in index_dict:
        row, col = index
        value = float(ui.parTW.item(row,col).text())
        i = name_list.index(index_dict[index])
        parameter_list[i][0] = value # update parameter value
        parameter_list[i][1] = False # clear the vary status first

    # update "vary" status in the selected cells
    selected_items = ui.parTW.selectionModel().selectedIndexes()
    selected_names = []
    for item in selected_items:
        selected_index = (item.row(), item.column())
        selected_name = index_dict[selected_index]
        i = name_list.index(selected_name)
        parameter_list[i][1] = True # update vary status
    
    return parameter_list
   
def name2index(name_list, reverse=False):
    ''' Create a mapping from table index to parameter name, or vice versa.
    If reverse=False, create a dict in the form of {name:position}, otherwise
    {position:name}.'''
    ndata = len([p for p in name_list if p.startswith('qoff')])
    nlayer = len([p for p in name_list if p.startswith('d')])
    
    index_dict = {}
    k = 0
    for row in range(ndata+nlayer+1):
        for col in range(4):
            if (row==0) & (col==0): continue
            elif (row>nlayer) & ((col==0)|(col==3)): continue
            else: 
                index_dict[name_list[k]] = (row,col)
                k = k + 1
    if reverse == False:
        return index_dict
    else:
        name_dict = {v:k for k,v in index_dict.iteritems()}
        return name_dict
  
def readData(rrf_files,sel_rows,fit_range,err_type=0):

    '''read multiple data set and cut them to fit range'''
    
    if sel_rows == []: return None
    if fit_range == None: fit_range = [0,1]
    # import pdb; pdb.set_trace();
    rrf = [np.loadtxt(rrf_files[r],comments='#') for r in sel_rows] 
    for i,d in enumerate(rrf):
        select = (d[:,0]>=fit_range[0])&(d[:,0]<=fit_range[1])
        rrf[i] = d[select]
    qz = tuple([a[:,0] for a in rrf]) # tuple: (qz1,qz2,...)
    y = tuple([a[:,1] for a in rrf]) # tuple: (data1,data2,...)
    if err_type==0: # tuple: (err1,err2,...)
        yerr = tuple([a[:,2] for a in rrf])
    elif err_type==1:
        yerr = tuple([np.sqrt(a[:,1]) for a in rrf])
    elif err_type==2:
        yerr = tuple([a[:,1] for a in rrf])
    else:
        yerr=tuple([np.ones(a[:,0].shape) for a in rrf])
       
    return (qz, y, yerr)
    
@jit(nopython=True)
def parratt(q, lamda, d, rho, beta):
    """
        Calculation of reflectivity by Parrat Recursion Formula without any roughness.
        Directly translated from a fortran subroutine 'parratt' in xr_ref.f90, which is
    attached as following:

    subroutine parratt(q,lambda,d,rho,beta,Rgen,Rgenr,M,N)
    !***************************************************************************
    !Calculation of reflectivity by Parratt Recursion Formula without any roughness
    !
    !M = No. of data points
    !N = No. of slabs
    !lambda = wavelength
    !d = list of thicknesses of each slab
    !rho=list of Electron densities of each slab
    !beta=list of Absorption coefficient in each slab
    !Rgen = generated reflectivtiy data
    !Rgenr= generated reflectance data
    !q = change in wave vector
    !***************************************************************************
    integer :: M,N
    double precision :: q(0:M), Rgen(0:M)
    double precision :: d(0:N+1), rho(0:N+1), beta(0:N+1), qc2(0:N+1)
    double precision :: lambda
    double complex :: X, fact1, fact2, r(0:N+1), k1, k2, fact,Rgenr(0:M)
    double precision, parameter :: re=2.817938e-5, pi=3.14159

    do j=0,N+1
       qc2(j)=16.0d0*pi*re*(rho(j)-rho(0))
    enddo

    do i = 0,M
       r(N+1)=dcmplx(0.0d0,0.0d0)
       do j=N,0,-1
          k1=cdsqrt(dcmplx(q(i)**2-qc2(j),-32.0d0*beta(j)*pi**2/lambda**2))
          k2=cdsqrt(dcmplx(q(i)**2-qc2(j+1),-32.0d0*beta(j+1)*pi**2/lambda**2))
          X=(k1-k2)/(k1+k2)
          fact1=dcmplx(dcos(dble(k2)*d(j+1)),dsin(dble(k2)*d(j+1)))
          fact2=dexp(-aimag(k2)*d(j+1))
          fact=fact1*fact2
          r(j)=(X+r(j+1)*fact)/(1.0+X*r(j+1)*fact)
       enddo
       Rgenr(i)=r(0)
       Rgen(i)=cdabs(r(0))**2
    enddo   
    end subroutine parratt

    """
    M = len(q) - 1
    N = len(d) - 2
    r = np.ones(N+2) * (0+0j) # this definition is compatible to numba
    qc2 = np.zeros(N+2)
    Rgen = np.zeros(M+1)
    Rgenr = np.ones(M+1) * (0+0j)

    for j in range(N+2):
        qc2[j] = 16.0 * np.pi * e_radius * (rho[j] - rho[0])
    

    for i in range(M+1):
        r[N+1] = 0+0j
        for j in range(N,-1,-1):
           k1 = np.sqrt(q[i]**2-qc2[j] - 32.0*beta[j]*np.pi**2/lamda**2*1j)
           k2 = np.sqrt(q[i]**2-qc2[j+1] - 32.0*beta[j+1]*np.pi**2/lamda**2*1j)
           X = (k1 - k2) / (k1 + k2)
           fac1 = m.cos(k2.real*d[j+1]) + m.sin(k2.real*d[j+1])*1j
           fac2 = m.exp(-k2.imag*d[j+1])
           fact = fac1 * fac2
           r[j] = (X + r[j+1]*fact) / (1.0+X*r[j+1]*fact)
        Rgenr[i] = r[0]
        Rgen[i] = abs(r[0])**2
    return (Rgen, Rgenr) 





################################################################################
################################################################################

if __name__ == '__main__':
    rho_t = 0.25913738441154344
    rho_b = 0.337
    itMu = 2.792661024598891e-09
    ibMu = 7.0179999999999995e-09
    qz = np.array([0.00879463, 0.00878271, 0.00877078, 0.00875884, 0.00874688,
       0.0087349 , 0.00872291, 0.0087109 , 0.00869887, 0.00868683,
       0.00867477, 0.00866269, 0.0086506 , 0.00863849, 0.00862636,
       0.00861422, 0.00860205, 0.00858987, 0.00857768, 0.00856546,
       0.00855323, 0.00854098, 0.00852872, 0.00851643, 0.00850413,
       0.00849181, 0.00847948, 0.00846712, 0.00845475, 0.00844236,
       0.00842995, 0.00841752, 0.00840507, 0.00839261, 0.00838012,
       0.00836762, 0.0083551 , 0.00834256, 0.00833   , 0.00831742,
       0.00830483, 0.00829221, 0.00827958, 0.00826692, 0.00825425,
       0.00824156, 0.00822884, 0.00821611, 0.00820336, 0.00819059,
       0.0081778 , 0.00816498, 0.00815215, 0.0081393 , 0.00812643,
       0.00811353, 0.00810062, 0.00808769, 0.00807473, 0.00806176,
       0.00804876, 0.00803574, 0.0080227 , 0.00800964, 0.00799656,
       0.00798346, 0.00797033, 0.00795719, 0.00794402, 0.00793083,
       0.00791762, 0.00790439, 0.00789113, 0.00787785, 0.00786455,
       0.00785123, 0.00783788, 0.00782452, 0.00781112, 0.00779771,
       0.00778427, 0.00777081, 0.00775733, 0.00774382, 0.00773029,
       0.00771673, 0.00770316, 0.00768955, 0.00767593, 0.00766227,
       0.0076486 , 0.0076349 , 0.00762118, 0.00760743, 0.00759365,
       0.00757985, 0.00756603, 0.00755218, 0.0075383 , 0.0075244 ,
       0.00751048, 0.00749652, 0.00748255, 0.00746854, 0.00745451,
       0.00744045, 0.00742637, 0.00741226, 0.00739812, 0.00738396,
       0.00736977, 0.00735555, 0.0073413 , 0.00732703, 0.00731272,
       0.00729839, 0.00728403, 0.00726965, 0.00725523, 0.00724079,
       0.00722631, 0.00721181, 0.00719728, 0.00718272, 0.00716813,
       0.00715351, 0.00713886, 0.00712418, 0.00710947, 0.00709473])
    ref = refCalFun([], [rho_t, rho_b], [itMu,ibMu], [3.0], qz)
    print(ref)
