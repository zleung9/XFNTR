import sys
# sys.path.append('/Users/zhuzi/work/data_analysis_20190514/')
import math as m
#import xr_ref as xr
import numpy as np
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
    k0=2*np.pi*float(energy)/12.3984 # wave vector
    theta=x/2/k0   # convert q to theta
   
    # total length of inner slabs plus 4 times roughness for both sides
    length = np.sum(d) + 4* (sigma[0]+sigma[-1])
    steps=int(length/slab) # each sliced box has thickness of ~ 0.25 A
    xsld=np.linspace(-4*sigma[0],np.sum(d)+4*sigma[-1],steps) # get the z-axis for sld

    sd=length/steps # thickness for each slab
    intrho=sldCalFun(d,rho,sigma,xsld) # rho for all the steps
    intmu=sldCalFun(d,mu,sigma,xsld) # mu for all the steps
    

    sdel=[]
    sbet=[]
    sdel.append(erad*2.0*np.pi/k0/k0*rho[0]) # delta for the top phase
    sbet.append(mu[0]/2/k0/1e8)        # beta for the top phase
    # add delta for the interface
    sdel=sdel+[intrho[i]*erad*2.0*np.pi/k0/k0 for i in range(len(intrho))] 
    sbet=sbet+[intmu[i]/2/k0/1e8 for i in range(len(intmu))] # add beta for the interface
    sdel.append(erad*2.0*np.pi/k0/k0*rho[-1])  # delta for the bottom phase
    sbet.append(mu[-1]/2/k0/1e8)    # beta for the bottom phase         
    
    d=slab*np.ones_like(sdel)
    lamda=2*np.pi/k0
    fdel=erad*2.0*np.pi/k0/k0
    sdelf=np.array(sdel)/fdel

    ref,refr=parratt(x,lamda,d,sdelf,sbet)
    frsnll,frsnl1=parratt(x,lamda,[0,1],[sdelf[0],sdelf[-1]],[sbet[0],sbet[-1]]) 
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
    r = np.zeros(N+2).astype(complex)
    qc2 = np.zeros(N+2)
    Rgen = np.zeros(M+1)
    Rgenr = np.zeros(M+1).astype(complex)

#    qc2 = 16 * np.pi * e_radius * (rho - rho[0]) 
    for j in range(N+2):
        qc2[j] = 16.0 * np.pi * e_radius * (rho[j] - rho[0])
    

    for i in range(M+1):
        for j in range(N+1)[::-1]:
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
    lamda = 0.61992
    d = np.array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])
    sdelf = np.array([0.25913738, 0.25913991, 0.25914098, 0.25914246, 0.25914451,
       0.25914731, 0.25915114, 0.2591563 , 0.25916323, 0.25917246,
       0.25918468, 0.25920073, 0.25922167, 0.25924878, 0.25928365,
       0.25932819, 0.25938466, 0.25945575, 0.25954462, 0.25965495,
       0.25979092, 0.25995733, 0.26015955, 0.26040355, 0.26069588,
       0.26104365, 0.26145444, 0.26193625, 0.26249736, 0.26314621,
       0.26389122, 0.26474062, 0.26570217, 0.266783  , 0.26798933,
       0.26932622, 0.27079734, 0.27240473, 0.27414862, 0.27602724,
       0.27803669, 0.28017091, 0.28242165, 0.2847785 , 0.28722903,
       0.28975897, 0.29235246, 0.29499231, 0.29766036, 0.3003379 ,
       0.30300595, 0.3056458 , 0.30823929, 0.31076923, 0.31321976,
       0.31557661, 0.31782735, 0.31996157, 0.32197102, 0.32384964,
       0.32559353, 0.32720092, 0.32867204, 0.33000893, 0.33121526,
       0.33229609, 0.33325764, 0.33410704, 0.33485205, 0.33550091,
       0.33606202, 0.33654382, 0.33695461, 0.33730238, 0.33759471,
       0.33783871, 0.33804093, 0.33820734, 0.33834332, 0.33845364,
       0.33854251, 0.33861361, 0.33867007, 0.33871461, 0.33874948,
       0.3387766 , 0.33879753, 0.33881358, 0.3388258 , 0.33883503,
       0.33884196, 0.33884712, 0.33885095, 0.33885375, 0.3388558 ,
       0.33885728, 0.33885835, 0.33886088])
    sbet = np.array([1.37766621e-18, 1.37778117e-18, 1.37782975e-18, 1.37789729e-18,
       1.37799053e-18, 1.37811835e-18, 1.37829232e-18, 1.37852743e-18,
       1.37884295e-18, 1.37926336e-18, 1.37981959e-18, 1.38055032e-18,
       1.38150350e-18, 1.38273809e-18, 1.38432587e-18, 1.38635345e-18,
       1.38892438e-18, 1.39216123e-18, 1.39620772e-18, 1.40123063e-18,
       1.40742155e-18, 1.41499818e-18, 1.42420521e-18, 1.43531447e-18,
       1.44862430e-18, 1.46445800e-18, 1.48316111e-18, 1.50509762e-18,
       1.53064483e-18, 1.56018692e-18, 1.59410732e-18, 1.63277989e-18,
       1.67655901e-18, 1.72576897e-18, 1.78069284e-18, 1.84156116e-18,
       1.90854091e-18, 1.98172519e-18, 2.06112403e-18, 2.14665680e-18,
       2.23814666e-18, 2.33531740e-18, 2.43779290e-18, 2.54509956e-18,
       2.65667159e-18, 2.77185929e-18, 2.88993998e-18, 3.01013148e-18,
       3.13160753e-18, 3.25351480e-18, 3.37499085e-18, 3.49518235e-18,
       3.61326304e-18, 3.72845074e-18, 3.84002277e-18, 3.94732943e-18,
       4.04980493e-18, 4.14697567e-18, 4.23846553e-18, 4.32399830e-18,
       4.40339713e-18, 4.47658141e-18, 4.54356117e-18, 4.60442948e-18,
       4.65935336e-18, 4.70856332e-18, 4.75234244e-18, 4.79101500e-18,
       4.82493541e-18, 4.85447750e-18, 4.88002470e-18, 4.90196122e-18,
       4.92066433e-18, 4.93649803e-18, 4.94980786e-18, 4.96091712e-18,
       4.97012415e-18, 4.97770078e-18, 4.98389170e-18, 4.98891461e-18,
       4.99296109e-18, 4.99619795e-18, 4.99876888e-18, 5.00079646e-18,
       5.00238424e-18, 5.00361883e-18, 5.00457201e-18, 5.00530274e-18,
       5.00585897e-18, 5.00627938e-18, 5.00659489e-18, 5.00683001e-18,
       5.00700398e-18, 5.00713179e-18, 5.00722504e-18, 5.00729258e-18,
       5.00734115e-18, 5.00745611e-18])
    x = np.array([0.00777308, 0.0077596 , 0.0077461 , 0.00773257, 0.00771902,
           0.00770544, 0.00769185, 0.00767822, 0.00766458, 0.00765091,
           0.00763721, 0.00762349, 0.00760974, 0.00759597, 0.00758218,
           0.00756836, 0.00755451, 0.00754064, 0.00752675, 0.00751282,
           0.00749888, 0.0074849 , 0.0074709 , 0.00745688, 0.00744282,
           0.00742874, 0.00741464, 0.00740051, 0.00738635, 0.00737216,
           0.00735794, 0.0073437 , 0.00732943, 0.00731513, 0.00730081,
           0.00728646, 0.00727207, 0.00725766, 0.00724322, 0.00722876,
           0.00721426, 0.00719973, 0.00718518, 0.00717059, 0.00715598,
           0.00714133, 0.00712666, 0.00711195, 0.00709722, 0.00708245,
           0.00706765, 0.00705282, 0.00703796, 0.00702307, 0.00700815,
           0.0069932 , 0.00697821, 0.00696319, 0.00694814, 0.00693306,
           0.00691794, 0.00690279, 0.00688761, 0.00687239, 0.00685714,
           0.00684185, 0.00682654, 0.00681118, 0.0067958 , 0.00678037,
           0.00676491, 0.00674942, 0.00673389, 0.00671833, 0.00670273,
           0.00668709, 0.00667141, 0.0066557 , 0.00663996, 0.00662417,
           0.00660835, 0.00659249, 0.00657659, 0.00656065, 0.00654467,
           0.00652866, 0.0065126 , 0.00649651, 0.00648037, 0.00646419,
           0.00644798, 0.00643172, 0.00641542, 0.00639908, 0.0063827 ,
           0.00636628, 0.00634981, 0.00633331, 0.00631675, 0.00630016,
           0.00628352, 0.00626684, 0.00625011, 0.00623334, 0.00621652,
           0.00619965, 0.00618274, 0.00616579, 0.00614879, 0.00613174,
           0.00611464, 0.00609749, 0.0060803 , 0.00606306, 0.00604576,
           0.00602842, 0.00601103, 0.00599359, 0.0059761 , 0.00595855,
           0.00594096, 0.00592331, 0.00590561, 0.00588785, 0.00587005,
           0.00585219, 0.00583427, 0.0058163 , 0.00579827, 0.00578019,
           0.00576205, 0.00574385, 0.00572559, 0.00570728, 0.0056889 ,
           0.00567047, 0.00565198, 0.00563343, 0.00561481, 0.00559614,
           0.0055774 , 0.00555859, 0.00553973, 0.0055208 , 0.0055018 ,
           0.00548274, 0.00546361, 0.00544442, 0.00542515, 0.00540582,
           0.00538642, 0.00536695, 0.00534741, 0.00532779, 0.00530811,
           0.00528835, 0.00526851, 0.00524861, 0.00522862, 0.00520856,
           0.00518842, 0.0051682 , 0.00514791, 0.00512753, 0.00510707,
           0.00508653, 0.00506591, 0.0050452 , 0.00502441, 0.00500353,
           0.00498256, 0.0049615 , 0.00494036, 0.00491912, 0.00489779,
           0.00487637, 0.00485485, 0.00483324, 0.00481153, 0.00478972,
           0.00476782, 0.00474581, 0.0047237 , 0.00470148, 0.00467916,
           0.00465673, 0.0046342 , 0.00461155, 0.00458879, 0.00456592,
           0.00454294, 0.00451983, 0.00449661, 0.00447327, 0.0044498 ,
           0.00442621, 0.0044025 , 0.00437865, 0.00435468, 0.00433057,
           0.00430633, 0.00428195, 0.00425743, 0.00423276, 0.00420796,
           0.00418301, 0.0041579 , 0.00413265, 0.00410724, 0.00408167,
           0.00405594, 0.00403004, 0.00400398, 0.00397775, 0.00395134,
           0.00392476, 0.00389799, 0.00387104, 0.0038439 , 0.00381657,
           0.00378904, 0.00376131, 0.00373337, 0.00370522, 0.00367686,
           0.00364827, 0.00361946, 0.00359042, 0.00356115, 0.00353162,
           0.00350186, 0.00347183, 0.00344154, 0.00341099, 0.00338016,
           0.00334904, 0.00331763, 0.00328593, 0.00325391, 0.00322158,
           0.00318891, 0.00315591, 0.00312256, 0.00308886, 0.00305477,
           0.00302031, 0.00298545, 0.00295017, 0.00291447, 0.00287832,
           0.00284172, 0.00280463, 0.00276705, 0.00272896, 0.00269032,
           0.00265112, 0.00261133, 0.00257093, 0.00252988, 0.00248816,
           0.00244572, 0.00240253, 0.00235855, 0.00231374, 0.00226804,
           0.0022214 , 0.00217376, 0.00212506, 0.00207521, 0.00202413,
           0.00197173, 0.0019179 , 0.00186251, 0.00180543, 0.00174648,
           0.00168547, 0.00162217, 0.00155629, 0.0014875 , 0.00141538,
           0.00133937, 0.00125878, 0.00117266, 0.0010797 , 0.00097795,
           0.00086429, 0.00073322, 0.00057292, 0.00034475])
    ref, refr = parratt(x,lamda,d,sdelf,sbet)
    print(ref)
    print('\n')
    print(refr)
