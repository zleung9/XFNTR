from collections import OrderedDict
import periodictable as pdtb
import lmfit as lm
import multiprocessing

r_e = pdtb.constants.electron_radius * 1e10  # classical electron radius, in A
N_A = pdtb.constants.avogadro_number  # Avogadro number, unitless
k_B = 1.38065e-23  # Boltzman constant, in J/K
p_igMu = (1/1.717)*1e-7    # absorption coefficient of 20keV beam by glass
p_thick = 3.5 * 1e7     # thickness of glass tray in A

# define an function that will be used in the program a lot.
absorb = lambda x: np.nan_to_num(np.exp(x))

import numpy as np
import scipy.stats as stat
import fit_ref as mfit
import flu_geometry_routines as gm



def penetrate(beta, delta, alpha, k0):
    alpha[alpha == np.inf] = 0
    alpha = alpha.astype(complex)
    beta_top, beta_bot = beta
    delta_top, delta_bot = delta
    alpha_c = np.sqrt(2 * (delta_bot - delta_top))
    trans = 4 * np.abs(alpha / (alpha + np.sqrt(alpha ** 2 - alpha_c ** 2))) ** 2
    penetration_coeff = 2 * k0 * np.imag(np.sqrt(alpha ** 2 - alpha_c ** 2 + beta_bot * 2j))
    return 1/penetration_coeff, trans

def update_flu_parameters(p, *args):

    assert type(p) is OrderedDict
    # *args have to be at least fitting parameters
    flu_par = args[0]

    # update the fitting parameter whatsoever
    try:
        p['hisc'] = flu_par['hisc'].value  # scale factor for the top phase, unitless.
        p['losc'] = flu_par['losc'].value  # scale factor for the bottom phase, unitless.
        p['bg'] = flu_par['bg'].value  # background intensity, unitless.
        p['tC'] = flu_par['upbk'].value  # ion concentration of the top phase, in M.
        p['bC'] = flu_par['lobk'].value  # ion concentration of the bottom phase, in M.
        p['sC'] = flu_par['surd'].value  # ion surface number density, in A^-2.
        p['qoff'] = flu_par['qoff'].value  # q off set for the data
        p['soff'] = flu_par['soff'].value * 1e7  # det range offset for the measurement
        p['l2off'] = flu_par['loff'].value * 1e7  # l2 offset for the measurement
        p['curv'] = flu_par['curv'].value * 1e10  # the curvature of the interface, in A.
        if p['curv'] == 0: p['curv'] = 10000 * 1e10
    except KeyError as e:
        print("Please check your parameter:{}".format(e))

    if len(args) == 1:
        return p

    # if the *args is tuple (flu_par, sys_par, flu_elements), do the following
    try:
        sys_par = args[1]
        flu_elements = args[2]
    except IndexError:
        print("update_flu_parameters takes 3 extra arguments!")

    # parameterize beam profile
    width = sys_par['width'] * 1e7  # width or FWHM of the beam, in A
    beam_profile = sys_par['beam']
    steps = 500
    if beam_profile == 'Uniform':
        beam_size = width
        weights = np.ones(steps + 1)
    elif beam_profile == 'Gaussian':
        stdev = width / 2.355  # FWHM of the beam, i.e. 2.355 sigma
        beam_size = 2 * (3 * stdev)  # keep the beam up to +/-3 standard deviation, or 99.73% of intensity.
        rays = np.linspace(-beam_size / 2, beam_size / 2, steps + 1)  # devide the beam into 500 rays
        weights = stat.norm(0, stdev).pdf(rays) * width  # weight normalized to the total intensity
    else:
        print("Error: please choose the right beam profile: uniform/gaussian")
        return None

    # unwrap system parameters
    E_inc = float(sys_par['E_inc'])  # energy of incidence, in KeV
    E_emit = float(sys_par['E_emt'])  # energy of  emission, in KeV
    mu_top_inc = float(sys_par['mu_top_inc'] / 1e8)  # abs. coef. of top phase for incidence, in 1/A
    mu_top_emit = float(sys_par['mu_top_emt'] / 1e8)  # abs. coef. of top phase for emission, in 1/A
    mu_bot_inc = float(sys_par['mu_bot_inc'] / 1e8)  # abs. coef. of bot phase for incidence, in 1/A
    mu_bot_emit = float(sys_par['mu_bot_emt'] / 1e8)  # abs. coef. of bot phase for emission, in 1/A
    rho_top = sys_par['rho_top']  # electron density of top pahse, in A^-3
    rho_bot = sys_par['rho_bot']  # electron density of top pahse, in A^-3

    det_len = sys_par['det_len'] * 1e7  # detector length, in A
    span = sys_par['span'] * 1e7

    # calculate abosrption parameters.
    k0 = 2 * np.pi * E_inc / 12.3984  # wave vector for incidence, in A^-1
    k1 = 2 * np.pi * E_emit / 12.3984  # wave vector for emission, in A^-1

    # construct elemental parameters
    vol_top, vol_bot = 0, 0  # vol_bot for ions for 1 L subphase
    ne_top, ne_bot = 0, 0  # total number of electrons from 1 L subphase
    bet_top_inc_ele, bet_top_emit_ele = 0, 0  # beta: for top phase at inc.& emit.
    bet_bot_inc_ele, bet_bot_emit_ele = 0, 0  # beta: for bot phase at inc.& emit.
    flupara_ele = {}  # setup dict for all elements in the subphase
    for i, e in enumerate(flu_elements):
        flupara_ele[i] = \
            [e[0], float(e[1]), float(e[2]),  # name, composition, ionic radius
             pdtb.elements.symbol(e[0]).number,
             pdtb.elements.symbol(e[0]).xray.scattering_factors(energy=sys_par['E_inc'])[1],
             pdtb.elements.symbol(e[0]).xray.scattering_factors(energy=sys_par['E_emt'])[1]]
    for i, p_ in flupara_ele.items():
        n_top_density = p['tC'] * p_[1] * N_A / 1e27  # atoms per A^3 in top phase
        n_bot_density = p['bC'] * p_[1] * N_A / 1e27  # atoms per A^3 in bot phase
        vol_top += n_top_density * 4 / 3 * np.pi * p_[2] ** 3
        vol_bot += n_bot_density * 4 / 3 * np.pi * p_[2] ** 3
        ne_top += n_top_density * p_[3]  # electrons per A^3
        ne_bot += n_bot_density * p_[3]  # electrons per A^3
        bet_top_inc_ele += n_top_density * 2 * np.pi * r_e * p_[4] / k0 ** 2
        bet_top_emit_ele += n_top_density * 2 * np.pi * r_e * p_[5] / k1 ** 2
        bet_bot_inc_ele += n_bot_density * 2 * np.pi * r_e * p_[4] / k0 ** 2
        bet_bot_emit_ele += n_bot_density * 2 * np.pi * r_e * p_[5] / k1 ** 2
    # absorption coefficient and electron density modified by solvent.
    rho_top = ne_top + (1 - vol_top) * rho_top
    rho_bot = ne_bot + (1 - vol_bot) * rho_bot

    # re-evaluate mu
    mu_top_inc = 2 * k0 * bet_top_inc_ele + (1 - vol_top) * mu_top_inc
    mu_top_emit = 2 * k1 * bet_top_emit_ele + (1 - vol_top) * mu_top_emit
    mu_bot_inc = 2 * k0 * bet_bot_inc_ele + (1 - vol_bot) * mu_bot_inc
    mu_bot_emit = 2 * k1 * bet_bot_emit_ele + (1 - vol_bot) * mu_bot_emit

    # calculate beta
    bet_top_inc = mu_top_inc / k0 / 2
    bet_top_emit = mu_top_emit / k1 / 2
    bet_bot_inc = mu_bot_inc / k0 / 2
    bet_bot_emit = mu_bot_emit / k1 / 2

    # calculate delta
    del_top_inc = 2 * np.pi * r_e * rho_top / k0 ** 2  # del=2*PI*re*rho/k^2, unitless #self.flutopdel
    del_top_emit = 2 * np.pi * r_e * rho_top / k1 ** 2
    del_bot_inc = 2 * np.pi * r_e * rho_bot / k0 ** 2
    del_bot_emit = 2 * np.pi * r_e * rho_bot / k1 ** 2

    p['k0'] = k0 # incident ray Energy, in KeV.
    p['k1'] = k1 # emission ray energy, in KeV.
    p['detR'] = det_len # detector range, in mm.
    p['wt'] = weights  # weights for the beam profile, the length of which is the total amount of steps.
    p['bmsz'] = beam_size  # the size of the beam for footprint calculation, in A.
    p['tRho'] = rho_top # electron density of top phase, in A^-3.
    p['bRho'] = rho_bot # electron density of bottom phase, in A^-3.
    p['itMu'] = mu_top_inc # mu for incident beam in top phase, in cm^-1
    p['etMu'] = mu_top_emit # mu for emitted beam in top phase, in cm^-1
    p['ibMu'] = mu_bot_inc # mu for incident beam in bottom phase, in cm^-1
    p['ebMu'] = mu_bot_emit # mu for emitted beam in bottom phase, in cm^-1
    p['itBt'] = bet_top_inc # beta for incident beam in top phase, in cm^-1
    p['etBt'] = bet_top_emit # beta for emitted beam in top phase, in cm^-1
    p['ibBt'] = bet_bot_inc # beta for incident beam in bottom phase, in cm^-1
    p['ebBt'] = bet_bot_emit # beta for emitted beam in bottom phase, in cm^-1
    p['itDt'] = del_top_inc # delta for incident beam in top phase, in cm^-1
    p['etDt'] = del_top_emit # delta for emitted beam in top phase, in cm^-1
    p['ibDt'] = del_bot_inc # delta for incident beam in bottom phase, in cm^-1
    p['ebDt'] = del_bot_emit # delta for emitted beam in bottom phase, in cm^-1
    p['span'] = span # the length of the sample cell, "the span", in A.

    return p

def fluCalFun_core(a0,sh,p):
    
    '''takes in flupara_fit, qz, return fluorescence data.
       Note that 'weights' contains the info of the steps for integration
       a0: the incident angle of X-ray beam, corrected with Qz_offset, in rad.
       sh: the sample height shift w.r.t. its norminal position, in A.
       p['wt']: weights for the beam profile, the length of which is the total amount of steps.
       p['k0']: wavevector for incident ray Energy, in KeV.
       p['detR']: detector range, in mm.
       p['hisc']: scale factor for the top phase, unitless.
       p['losc']: scale factor for the bottom pahse, unitless.
       p['bg']: background intensity, unitless.
       p['k1']: wavevector for emission ray energy, in KeV.
       p['tC']: ion concentration of the top phase, in M.
       p['bC']: ion concentration of the bottom phase, in M.
       p['sC']: ioin surface number density, in A^-2.
       p['qoff']: q off set for the data
       p['soff']: sample height offset (effoct of detector offset included) for the measurement
       p['l2off']: l2 offset for the measurement
       p['tRho']: electron density of top phase, in A^-3.
       p['bRho']: electron density of bottom phase, in A^-3.
       p['itMu']: mu for incident beam in top pahse, in A^-1
       p['etMu']: mu for emitted beam in top phass, in A^-1
       p['ibMu']: mu for incident beam in bottom pahse, in A^-1
       p['ebMu']: mu for emitted beam in bottom phass, in A^-1
       p['itBt']: beta for incident beam in top pahse, in A^-1
       p['etBt']: beta for emitted beam in top phass, in A^-1
       p['ibBt']: beta for incident beam in bottom pahse, in A^-1
       p['ebBt']: beta for emitted beam in bottom phass, in A^-1
       p['itDt']: delta for incident beam in top pahse, in A^-1
       p['etDt']: delta for emitted beam in top phass, in A^-1
       p['ibDt']: delta for incident beam in bottom pahse, in cm^-1
       p['ebDt']: delta for emitted beam in bottom phass, in cm^-1
       p['span']: the length of the sample cell, "the span", in A.
       p['curv']: the curvature of the interface, in A.
       p['bmsz']: the size of the beam for footprint calculation, in A.
       '''


    steps = len(p['wt']) - 1
    fprint = p['bmsz'] / np.sin(a0) # foortprint of the beam on the interface.
    stepsize = fprint / steps
    center = - sh / a0

    # initialize fluorescence data, rows: total, aqueous, organic, interface


    # get the position of single ray hitting the surface
    x0 = np.linspace(center-fprint/2, center+fprint/2, steps+1)
    surface = np.array([gm.hit_surface([xx, 0], -a0, p['curv'], p['span']) for xx in x0])
    hit = np.isfinite(surface[:,0]) * np.isfinite(surface[:,1])  # rays that hit the liquid-liquid interface
    block = np.isfinite(surface[:,0]) * np.isinf(surface[:,1])  # rays that are blocked by the tray
    x_s = surface[:, 0][hit]  # x' for rays that hit on the interface
    z_s = surface[:, 1][hit] # z' for rays that hit on the interface
    wt_s = p['wt'][hit]  # weight for rays that hit on the interface

    # (x,z) and other surface geometry for points where beam hit at the interface.
    theta = -x_s / p['curv']  # incident angle
    a_new = a0 + theta  # actual incident angle w.r. to the surface
    # a1 = a0 + 2 * theta
    a1 = a0
    x_inc = x_s + z_s / a0  # x position where the inc. xray passes z=0 line
    x_ref = x_s - z_s / a1  # x position where the ref. xray passes z=0 line.

    mu_eff_inc = p['etMu'] + p['itMu'] / a0  # eff.abs.depth for incident beam in oil phase
    mu_eff_ref = p['etMu'] - p['itMu'] / a1  # eff.abs.depth for reflected beam in water phase

    # z coordinate of the intersection of ray with following:
    z_inc_l = (x_inc + p['detR'] / 2) * a0  # incidence with left det. boundary: x=-l/2
    z_inc_r = (x_inc - p['detR'] / 2) * a0  # incidence with right det. boundary: x=l/2
    z_ref_l = -(x_ref + p['detR'] / 2) * a1  # reflection with left det. boundary: x=-l/2
    z_ref_r = -(x_ref - p['detR'] / 2) * a1  # reflection with right det. boundary: x=l/2

    # two regions: region3: [-h/2a0,-l/2] & region 2: [-l/2,l/2]
    x_region = [(x_s <= -p['detR'] / 2), (x_s > -p['detR'] / 2) * (x_s < p['detR'] / 2)]

    ################### for redgion 1: region x>= l/2  ########################
    x0_region1 = x0[surface[:, 0] > p['detR'] / 2]  # choose x0 with x'>l/2
    wt_region1 = p['wt'][surface[:, 0] > p['detR'] / 2]  # choose weight with x'>l/2
    upper_bulk1 = wt_region1 * \
                  absorb(-x0_region1 * p['itMu']) / mu_eff_inc * \
                  (absorb((x0_region1 + p['detR'] / 2) * a0 * mu_eff_inc) -
                   absorb((x0_region1 - p['detR'] / 2) * a0 * mu_eff_inc))

    # define the initial intensity to be just background value
    flu = np.array([p['bg']] * 6)

    # if beam miss the surface entirely, do the following:
    if len(x_s) == 0:  # the entire beam miss the interface, only incidence in upper phase.
        # sh_offset_factor = absorb(-mu_top_emit * center[i] * a0)
        usum_inc = stepsize * np.sum(upper_bulk1)
        flu[3] += p['hisc'] * usum_inc * N_A * p['tC'] / 1e27  # oil phase incidence only + background
        flu[0] = flu[3]  # total intensity only contains oil phase
        return flu

    ref = mfit.refCalFun([], [p['tRho'], p['bRho']], [p['itMu'], p['ibMu']], [3.0], 2 * p['k0'] * a_new)
    p_depth, trans = penetrate((p['itBt'],p['ibBt']), (p['itDt'],p['ibDt']), a_new, p['k0'])
    p_depth_eff = 1 / (p['ebMu'] + a_new/a0 / p_depth)

    ################### for region -l/2 < x < l/2  #################
    lower_bulk2 = x_region[1] * wt_s * absorb(-x_s * p['itMu'] - z_s / p_depth) * trans * p_depth * \
                  (absorb(z_s / p_depth_eff) - absorb(z_inc_r / p_depth_eff))
    interface = x_region[1] * wt_s * trans * absorb(-p['itMu'] * x_s)
    upper_bulk2_inc = x_region[1] * wt_s * \
                      (absorb(-x_inc * p['itMu']) / mu_eff_inc * (
                              absorb(z_inc_l * mu_eff_inc) - absorb(z_s * mu_eff_inc)))
    upper_bulk2_inc[np.isnan(upper_bulk2_inc)] = 0  # if there is nan, set to 0
    upper_bulk2_ref = x_region[1] * wt_s * \
                      (absorb(-x_ref * p['itMu']) / mu_eff_ref * ref * (
                              absorb(z_ref_r * mu_eff_ref) - absorb(z_s * mu_eff_ref)))
    upper_bulk2_ref[np.isnan(upper_bulk2_ref)] = 0  # if there is nan, set to 0

    ###################### for region x<=-l/2 ########################
    lower_bulk3 = x_region[0] * wt_s * absorb(-x_s * p['itMu'] - z_s / p_depth) * trans * p_depth_eff * \
                  (absorb(z_inc_l / p_depth_eff) - absorb(z_inc_r / p_depth_eff))
    upper_bulk3 = x_region[0] * wt_s * absorb(-x_ref * p['itMu']) / mu_eff_ref * ref * \
                  (absorb(mu_eff_ref * z_ref_r) - absorb(mu_eff_ref * z_ref_l))
    # if there are rays that are blocked by tray, their intensity is still significent
    if np.sum(block) > 0:
        edge = None
    # combine the two regions and integrate along x direction by performing np.sum.
    bsum = stepsize * np.sum(lower_bulk3 + lower_bulk2)
    ssum = stepsize * np.sum(interface)
    usum_inc = stepsize * (np.sum(upper_bulk1) + np.sum(upper_bulk2_inc))
    usum_ref = stepsize * (np.sum(upper_bulk3) + np.sum(upper_bulk2_ref))


    # vectorized integration method is proved to reduce the computation time by a factor of 5 to 10.
    int_bulk = p['losc'] * bsum * N_A * p['bC'] / 1e27
    int_upbk_inc = p['hisc'] * usum_inc * N_A * p['tC'] / 1e27  # metal ions in the upper phase.
    int_upbk_ref = p['hisc'] * usum_ref * N_A * p['tC'] / 1e27  # metal ions in the upper phase.
    int_sur = p['losc'] * ssum * p['sC']

    flu += np.array([int_bulk+int_sur+int_upbk_inc+int_upbk_ref,  # 2. total fluorescence + background
                     int_bulk,  # 3. lower bulk
                     int_sur,  # 4. interface
                     int_upbk_inc+int_upbk_ref, # 5. upper bulk
                     int_upbk_inc,  # 6. upper bulk incidence
                     int_upbk_ref])  # 7. upper bulk reflection
    return flu

def flu2min(pars, x, p, data=None, eps=None): # residuel for flu fitting

    sh, qz = x
    p = update_flu_parameters(p, (pars))

    alpha = (qz + p['qoff']) / p['k0'] / 2 # include the qz offset.
    a0006 = 0.006 / p['k0'] / 2  # incident angle for qz=0.006

    # initialize fluorescence data 3-D matrix
    flu = np.zeros((len(sh), len(qz), 8))
    try:
        for i, ds in enumerate(sh):
            for j, a0 in enumerate(alpha):
                dsh = -p['l2off'] * (a0 - a0006) + p['soff'] + ds*1e7
                flu[i, j, 0] = ds
                flu[i, j, 1] = qz[j]
                flu[i, j, 2:] = fluCalFun_core(a0, dsh, p)
    except KeyError as e:
        print("Please check parameter: {}".format(e))
    if data is None:
        return flu
    if eps is None:
        return (flu[:,:,2].flatten() - data)

    return (flu[:,:,2].flatten() - data) / eps

def fluErrorFitSingleCore2(i, value_list, sh, qz, pname, flu_par, flucal_par, data_to_fit):
    flu_par[pname].value = value_list[i]  # change value of the chosen parameter
    fluerr_result = lm.minimize(fl.flu2min, flu_par,
                                args=((sh, qz), flucal_par),
                                kws={'data': data_to_fit[:, 1], 'eps': data_to_fit[:, 2]})
    return ([i, value_list[i], fluerr_result.nfree, fluerr_result.redchi])

def multiCore(func, iterable):
    pool = multiprocessing.Pool()
    result = pool.map(func, range(len(iterable)))
    pool.close()
    pool.join()
    return result

if __name__ == '__main__':
    p_depth, trans = penetrate((1.35e-10, 5.0e-10), (4.466e-7, 5.84e-7), np.array([3.14]), 1.1355)
    print(p_depth)
    print(trans)