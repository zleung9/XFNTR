#import xr_ref as xr
import sys
import math as m
import numpy as np
from numba import jit
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
    sbet = np.array(sbet)
    
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

#    qc2 = 16 * np.pi * e_radius * (rho - rho[0]) 
    for j in range(N+2):
        qc2[j] = 16.0 * np.pi * e_radius * (rho[j] - rho[0])
    

    for i in range(M+1):
        for j in range(N+1):
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
    d = np.array([])
    d = np.array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])
    sdelf =np.array([0.25913738, 0.25913972, 0.25914071, 0.25914209, 0.25914398,
       0.25914658, 0.25915013, 0.25915491, 0.25916133, 0.25916988,
       0.2591812 , 0.25919607, 0.25921547, 0.25924059, 0.2592729 ,
       0.25931416, 0.25936648, 0.25943234, 0.25951469, 0.2596169 ,
       0.25974288, 0.25989705, 0.26008441, 0.26031047, 0.26058131,
       0.26090351, 0.2612841 , 0.26173049, 0.26225035, 0.2628515 ,
       0.26354175, 0.26432869, 0.26521956, 0.26622093, 0.26733858,
       0.26857719, 0.26994016, 0.27142938, 0.27304507, 0.27478558,
       0.27664731, 0.27862463, 0.28070991, 0.28289349, 0.28516387,
       0.28750783, 0.28991065, 0.29235643, 0.29482835, 0.29730904,
       0.29978095, 0.30222673, 0.30462956, 0.30697351, 0.30924389,
       0.31142748, 0.31351275, 0.31549008, 0.31735181, 0.31909231,
       0.320708  , 0.32219723, 0.3235602 , 0.32479881, 0.32591645,
       0.32691783, 0.32780869, 0.32859564, 0.32928589, 0.32988704,
       0.3304069 , 0.33085328, 0.33123387, 0.33155607, 0.33182692,
       0.33205298, 0.33224033, 0.33239451, 0.33252049, 0.3326227 ,
       0.33270504, 0.33277091, 0.33282322, 0.33286448, 0.33289679,
       0.33292191, 0.33294131, 0.33295618, 0.3329675 , 0.33297605,
       0.33298248, 0.33298726, 0.3329908 , 0.3329934 , 0.3329953 ,
       0.33299667, 0.33299766, 0.333     ])
    sbet = np.array([1.37766621e-18, 1.37773223e-18, 1.37776013e-18, 1.37779891e-18,
       1.37785246e-18, 1.37792586e-18, 1.37802576e-18, 1.37816078e-18,
       1.37834196e-18, 1.37858339e-18, 1.37890280e-18, 1.37932243e-18,
       1.37986980e-18, 1.38057877e-18, 1.38149056e-18, 1.38265492e-18,
       1.38413129e-18, 1.38599007e-18, 1.38831379e-18, 1.39119823e-18,
       1.39475340e-18, 1.39910433e-18, 1.40439153e-18, 1.41077109e-18,
       1.41841435e-18, 1.42750695e-18, 1.43824733e-18, 1.45084452e-18,
       1.46551516e-18, 1.48247990e-18, 1.50195891e-18, 1.52416688e-18,
       1.54930732e-18, 1.57756645e-18, 1.60910684e-18, 1.64406086e-18,
       1.68252441e-18, 1.72455096e-18, 1.77014625e-18, 1.81926399e-18,
       1.87180263e-18, 1.92760354e-18, 1.98645076e-18, 2.04807229e-18,
       2.11214325e-18, 2.17829052e-18, 2.24609912e-18, 2.31511986e-18,
       2.38487826e-18, 2.45488429e-18, 2.52464269e-18, 2.59366343e-18,
       2.66147202e-18, 2.72761929e-18, 2.79169025e-18, 2.85331179e-18,
       2.91215900e-18, 2.96795992e-18, 3.02049856e-18, 3.06961630e-18,
       3.11521159e-18, 3.15723813e-18, 3.19570168e-18, 3.23065571e-18,
       3.26219609e-18, 3.29045523e-18, 3.31559567e-18, 3.33780364e-18,
       3.35728265e-18, 3.37424738e-18, 3.38891803e-18, 3.40151521e-18,
       3.41225559e-18, 3.42134819e-18, 3.42899145e-18, 3.43537101e-18,
       3.44065821e-18, 3.44500914e-18, 3.44856431e-18, 3.45144875e-18,
       3.45377247e-18, 3.45563126e-18, 3.45710763e-18, 3.45827198e-18,
       3.45918377e-18, 3.45989274e-18, 3.46044012e-18, 3.46085974e-18,
       3.46117916e-18, 3.46142058e-18, 3.46160177e-18, 3.46173679e-18,
       3.46183669e-18, 3.46191009e-18, 3.46196363e-18, 3.46200242e-18,
       3.46203031e-18, 3.46209633e-18])
    x = np.array([0.00923793, 0.00921903, 0.0092001 , 0.00918112, 0.00916211,
       0.00914306, 0.00912396, 0.00910483, 0.00908566, 0.00906644,
       0.00904719, 0.00902789, 0.00900855, 0.00898917, 0.00896975,
       0.00895029, 0.00893079, 0.00891124, 0.00889165, 0.00887201,
       0.00885233, 0.00883261, 0.00881285, 0.00879304, 0.00877318,
       0.00875328, 0.00873334, 0.00871335, 0.00869331, 0.00867323,
       0.0086531 , 0.00863292, 0.00861269, 0.00859242, 0.0085721 ,
       0.00855174, 0.00853132, 0.00851085, 0.00849034, 0.00846977,
       0.00844916, 0.0084285 , 0.00840778, 0.00838701, 0.00836619,
       0.00834532, 0.0083244 , 0.00830343, 0.0082824 , 0.00826131,
       0.00824018, 0.00821899, 0.00819774, 0.00817644, 0.00815509,
       0.00813367, 0.00811221, 0.00809068, 0.0080691 , 0.00804746,
       0.00802576, 0.008004  , 0.00798218, 0.00796031, 0.00793837,
       0.00791637, 0.00789431, 0.00787219, 0.00785001, 0.00782776,
       0.00780545, 0.00778307, 0.00776064, 0.00773813, 0.00771556,
       0.00769293, 0.00767023, 0.00764746, 0.00762462, 0.00760171,
       0.00757874, 0.00755569, 0.00753258, 0.00750939, 0.00748613,
       0.0074628 , 0.0074394 , 0.00741592, 0.00739237, 0.00736874,
       0.00734503, 0.00732125, 0.00729739, 0.00727346, 0.00724944,
       0.00722535, 0.00720117, 0.00717691, 0.00715257, 0.00712815,
       0.00710364, 0.00707905, 0.00705437, 0.00702961, 0.00700476,
       0.00697982, 0.00695479, 0.00692967, 0.00690446, 0.00687916,
       0.00685376, 0.00682827, 0.00680268, 0.006777  , 0.00675121,
       0.00672533, 0.00669935, 0.00667327, 0.00664709, 0.0066208 ,
       0.00659441, 0.00656791, 0.00654131, 0.00651459, 0.00648777,
       0.00646083, 0.00643378, 0.00640662, 0.00637934, 0.00635195,
       0.00632444, 0.0062968 , 0.00626905, 0.00624117, 0.00621316,
       0.00618503, 0.00615677, 0.00612838, 0.00609986, 0.0060712 ,
       0.00604241, 0.00601348, 0.00598441, 0.0059552 , 0.00592584,
       0.00589634, 0.00586669, 0.00583689, 0.00580694, 0.00577683,
       0.00574656, 0.00571613, 0.00568554, 0.00565479, 0.00562387,
       0.00559277, 0.0055615 , 0.00553006, 0.00549843, 0.00546662,
       0.00543463, 0.00540245, 0.00537007, 0.0053375 , 0.00530473,
       0.00527175, 0.00523856, 0.00520517, 0.00517156, 0.00513773,
       0.00510367, 0.00506939, 0.00503487, 0.00500011, 0.00496511,
       0.00492987, 0.00489436, 0.0048586 , 0.00482258, 0.00478628,
       0.0047497 , 0.00471285, 0.0046757 , 0.00463825, 0.0046005 ,
       0.00456244, 0.00452405, 0.00448534, 0.00444629, 0.00440689,
       0.00436714, 0.00432703, 0.00428654, 0.00424566, 0.00420438,
       0.0041627 , 0.00412059, 0.00407805, 0.00403506, 0.00399161,
       0.00394768, 0.00390326, 0.00385832, 0.00381285, 0.00376684,
       0.00372026, 0.00367308, 0.00362529, 0.00357687, 0.00352778,
       0.00347799, 0.00342748, 0.00337622, 0.00332417, 0.00327129,
       0.00321754, 0.00316287, 0.00310725, 0.00305061, 0.0029929 ,
       0.00293405, 0.002874  , 0.00281267, 0.00274997, 0.00268581,
       0.00262007, 0.00255265, 0.00248339, 0.00241215, 0.00233874,
       0.00226294, 0.00218452, 0.00210318, 0.00201856, 0.00193023,
       0.00183767, 0.00174019, 0.00163691, 0.00152666, 0.00140781,
       0.00127795, 0.0011333 , 0.00096726, 0.00076604])
    ref, refr = parratt(x,lamda,d,sdelf,sbet)
    print(ref)
    print('\n')
    print(refr)
