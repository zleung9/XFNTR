import sys
# sys.path.append('/Users/zhuzi/work/data_analysis_20190514/')
import xr_ref as xr
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

def refCalFun(d,rho,mu,sigma,syspara,x,rrf=True):

    if sigma[0] <= 0: sigma[0] = 1e-5 # eliminate zero
    sigma = sigma[0] * np.ones(len(sigma)) # fixed sigma

    qoff = syspara[0]
    yscale = syspara[1]
    qres = syspara[2]
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

    ref,refr=xr.parratt(x+qoff,lamda,d,sdelf,sbet)
    frsnll,frsnl1=xr.parratt(x,lamda,[0,1],[sdelf[0],sdelf[-1]],[sbet[0],sbet[-1]]) 
    if rrf == True:
        return yscale*ref/frsnll
    else:
        return yscale*ref

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
                           [qoff[i]]+[1,0],
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
                           [qoff[i]] + [1, 0],
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
        print iteration,redchisq

def initParameters(name_list,para_list):
    
    if len(name_list)!=len(para_list):
        print " Does name_list and value_list match, please check "
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
    
################################################################################
################################################################################

if __name__ == '__main__':
    
    # sys.path.append('/Users/zhuzi/work/data_analysis_20180704/')
    
    true_params = Parameters()
    true_params.add_many( 
    
            # (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
            # different parameters
            ('rho_b1', 0.33, False, 0, None, None, None),
            ('mu_b1',      0, False, 0, None, None, None),
            ('qoff1',-5.05e-4, True, None, None, None, None), # To be fitted
            ('rho_b2', 0.35, False, 0, None, None, None), 
            ('mu_b2',      0, False, 0, None, None, None),
            ('qoff2', 5.05e-4, True, None, None, None, None), # To be fitted
            ('rho_b3', 0.37, False, 0, None, None, None),
            ('mu_b3',      0, False, 0, None, None, None),
            ('qoff3',-3.05e-4, True, None, None, None, None), # To be fitted
            ('rho_b4', 0.39, False, 0, None, None, None),
            ('mu_b4',      0, False, 0, None, None, None),
            ('qoff4', 7.05e-4, True, None, None, None, None), # To be fitted
            ('rho_b5', 0.43, False, 0, None, None, None),
            ('mu_b5',      0, False, 0, None, None, None),
            ('qoff5', 3.05e-4, True, None, None, None, None), # To be fitted
        
            # shared parameters
            ('rho_t', 0.2592, False, None, None, None, None),
            ('mu_t',       0, False, None, None, None, None),
            ('d1',      24, True, 0, None, None, None), # To be fitted
            ('rho1', 0.355, True, 0, None, None, None), # To be fitted
            ('mu1',      0, False, 0, None, None, None),
            ('sigma0',    2, True, 0, None, None, None), # To be fitted
            ('y_scale',   1, False, 0, None, None, None),
            ('q_res',     0, False, 0, None, None, None),
    )
    
    file_dir = "/Users/zhuzi/Documents/2018_b_summer/research/agent/"

    file_rrf = ["simulated_reflectivity_1_rrf.txt",
                "simulated_reflectivity_2_rrf.txt",
                "simulated_reflectivity_3_rrf.txt",
                "simulated_reflectivity_4_rrf.txt",
                "simulated_reflectivity_5_rrf.txt"]

    file_fit = ["simulated_reflectivity_1_fit.txt",
                "simulated_reflectivity_2_fit.txt",
                "simulated_reflectivity_3_fit.txt",
                "simulated_reflectivity_4_fit.txt",
                "simulated_reflectivity_5_fit.txt"]

    rrf = [np.loadtxt(file_dir+f) for f in file_rrf]

    qz = tuple([a[:,0] for a in rrf])
    data = tuple([a[:,1] for a in rrf])
    err = tuple([a[:,2] for a in rrf])
    source = ref2min(true_params,qz,data,err,fit=False)

    ############################################################################
    params = Parameters()
    params.add_many( # (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)

                # different parameters
                ('rho_b1', 0.333, False, 0, None, None, None),
                ('mu_b1',      0, False, 0, None, None, None),
                ('qoff1',-5.05e-4, True, None, None, None, None), # To be fitted
                ('rho_b2', 0.35, False, 0, None, None, None),
                ('mu_b2',      0, False, 0, None, None, None),
                ('qoff2', 5.05e-4, True, None, None, None, None), # To be fitted
                ('rho_b3', 0.37, False, 0, None, None, None),
                ('mu_b3',      0, False, 0, None, None, None),
                ('qoff3',-3.05e-4, True, None, None, None, None), # To be fitted
                ('rho_b4', 0.39, False, 0, None, None, None),
                ('mu_b4',      0, False, 0, None, None, None),
                ('qoff4', 7.05e-4, True, None, None, None, None), # To be fitted
                ('rho_b5', 0.433, False, 0, None, None, None),
                ('mu_b5',      0, False, 0, None, None, None),
                ('qoff5', 3.05e-4, True, None, None, None, None), # To be fitted

                # shared parameters
                ('rho_t', 0.2592, False, None, None, None, None),
                ('mu_t',       0, False, None, None, None, None),
                ('d1',      24, True, 0, None, None, None), # To be fitted
                ('rho1', 0.355, True, 0, None, None, None), # To be fitted
                ('mu1',      0, False, 0, None, None, None),
                ('sigma0',    2, True, 0, None, None, None), # To be fitted
                ('y_scale',   1, False, 0, None, None, None),
                ('q_res',     0, False, 0, None, None, None),
    )

    guess = ref2min(params,qz,data,err,fit=False)
    out = minimize(ref2min, params, 
                   args=(qz,data,err),
                   kws={'fit':True})

    # for p in params: print params[p]
    print '\n'
    print out.redchi
    report_fit(out)
    ############################################################################
    fit = ref2min(out.params,qz,data,err,fit=False)

    fig,ax = plt.subplots(figsize=(10,6))
    for i in range(len(qz)):
        if i in [0,1,2,3,4]:
            ax.errorbar(qz[i],data[i],yerr=err[i],
                        ls='', marker='o',color='r',
                        label='data'+str(i+1))
            ax.plot(qz[i],fit[i],color='g',
                    label='data source')
    ax.legend(loc='upper right')
    ax.grid(1)
    plt.show()