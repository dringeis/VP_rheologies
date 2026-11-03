#!/usr/bin python
import numpy as np
# from smoothabs import smoothmin_good
from matplotlib.patches import Ellipse
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
from matplotlib.colors import SymLogNorm
from math import copysign

# import smoothmin functions
from smooth_min import *

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

from VP_rheology_settings import *

def create_data(random=True,i=1e-6,j=0,plot=False,sym=True,s=201):
    '''
    Creating fake random data
    note that this not the correct way to compute the strain rates on a C grid
    but it is irrelevant for the goal of this code
    '''
    if random :

        # change the seed to vary the results
        np.random.seed(42)
        u=np.random.random((s,s))*i+j
        np.random.seed(24)
        v=np.random.random((s,s))*i+j

        #  For the mu(I) rheology
        np.random.seed(4)
        phi=np.random.random((s-1,s-1))
        np.random.seed(2)
        h=np.random.random((s-1,s-1))

        # spatial grid spacing
        dxy=1000

        # deformations
        dudx=(u[1:,1:]-u[:-1,1:])/dxy
        dvdx=(v[1:,1:]-v[:-1,1:])/dxy
        dudy=(u[1:,1:]-u[1:,:-1])/dxy
        dvdy=(v[1:,1:]-v[1:,:-1])/dxy

        e11=dudx
        e22=dvdy
        if sym:
            e12=0.5*(dudy+dvdx)
            e21=e12
        else:
            e12=dudy
            e21=dvdx

    else:
        eg=np.mgrid[-i:i:200j,i:-i:200j]
        e11 = eg[1,:,:]#*0.0
        e22 = eg[0,:,:]#*-0.5
        e12 = (1.0*e11+1.0*e22)*0.0
        e21 = (1.0*e11+1.0*e22)*0.0

        #  For the mu(I) rheology
        phi = np.ones((200,200))
        h = np.ones((200,200))

    if plot :
        plt.figure('e11')
        plt.pcolormesh(e11)
        plt.colorbar()
        plt.axis('equal')

        plt.figure('e22')
        plt.pcolormesh(e22)
        plt.colorbar()
        plt.axis('equal')

        plt.figure('e12')
        plt.pcolormesh(e12)
        plt.colorbar()
        plt.axis('equal')

        if not sym :
            plt.figure('e21')
            plt.pcolormesh(e21)
            plt.colorbar()
            plt.axis('equal')

        plt.show()



    return {'e11':e11, 'e22':e22, 'e12':e12, 'e21':e21, 'h':h, 'phi':phi}

def comp_sim_sr(data):
    '''
    Computing simpler variables
    '''

    data['ep'] = data['e11']+data['e22']
    data['em'] = data['e11']-data['e22']
    data['e12s'] = 0.5*(data['e12']+data['e21'])

    ### Computing deformations in principal stress
    # eTmp=np.sqrt(em**2+4*np.abs(e12*e21))
    data['eTmp'] = np.sqrt(np.abs(data['em']**2+4*data['e12']*data['e21']))
    # eTmp=np.sqrt(em**2+4*e12*e21)

    data['e1'] = 0.5*(data['ep']+data['eTmp'])
    data['e2'] = 0.5*(data['ep']-data['eTmp'])

    ### Computing in stress invariant
    data['eI'] = data['e1']+data['e2']
    data['eII'] = data['e1']-data['e2']

    data['eI'][data['eI']==0] = 1e-10
    data['eII'][data['eII']==0] = 1e-10

    ### Computing deformations in principal stress
    data['eTmps'] = np.sqrt(data['em']**2+4*data['e12s']**2)
    data['e1s'] = 0.5*(data['ep']+data['eTmp'])
    data['e2s'] = 0.5*(data['ep']-data['eTmp'])

    ### Computing in stress invariant
    data['eIs'] = 0.5*(data['e1']+data['e2'])
    data['eIIs'] = 0.5*(data['e1']-data['e2'])

    # return ep, em, e1, e2, eI, eII, e12s, e1s, e2s, eIs, eIIs
    return None


##########################
# RHEOLOGIES (VISCOSITIES)
##########################

def compute_visc(data={},rheos={}):

    for rheo_n in rheos:
        data['rheos'].append(rheo_n)
        data[rheo_n] = {}

        if rheos[rheo_n]['rheo_t'] == 'ell' :
            ellip(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'ellt' :
            ellip_test(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'ell_rot' :
            ellip(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n, rot=True)
        elif rheos[rheo_n]['rheo_t'] == 'mce' :
            mce(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mceG' :
            mceG(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mcs' :
            mcs(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mcpl' :
            mcpl(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mctd' :
            mctd(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'pl' :
            pl(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'td' :
            td(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'epl' :
            epl(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'etd' :
            etd(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'muID' :
            muID(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)

    return None

def ellip(data={}, rheo={}, rheo_n = '', rot=False):
    '''
    ELLIPTICAL YIELD CURVE
    '''
    print('Computing Ellipse rheology ; name:',rheo_n)


    # load data
    ep = data['ep']
    em = data['em']
    if rot :
        e12 = data['e12']
        e21 = data['e21']
    else:
        e12 = data['e12s']
        e21 = data['e12s']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'efr' in rheo:
        efr = rheo['efr']
    elif e != 2:
        efr = e
        rheo['efr'] = efr
    else:
        efr = e_d
        rheo['efr'] = efr

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    recip_e2 = 1./(e*e)
    recip_efr2 = 1./(efr*efr)
    e2_recip_efr4 = e**2/(efr**4)

    ### Computing Delta
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*np.abs(e12*e21))
    # deltaCsq=ep*ep+e2_recip_efr4*(np.abs(em*em+4.0*e12*e21))
    deltaCsq = ep**2 + e2_recip_efr4 * (em**2 + 4.0*e12**2)
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*e12*e21)

    deltaC = np.sqrt(deltaCsq)
    ## with a sqrt function
    if SEAICE_DELTA_SMOOTHREG :
        deltaCreg = np.sqrt(deltaCsq + deltaMinSq)
    ## with a sharp min function
    else:
        deltaCreg = np.sqrt(np.maximum(deltaCsq,deltaMinSq))

    ### Computing Zeta
    ## with a smooth tanh function
    if SEAICE_ZETA_SMOOTHREG:
        argTmp = np.exp(-1./(deltaCreg*SEAICE_zetaMaxFac))
        zeta=ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+kt)
    ## with a sharp min function
    else:
        zeta=press0*(1.+kt)/(2.*deltaCreg)
        zeta=np.minimum(zeta,ZMAX)
    zeta=np.maximum(zeta,ZMIN)

    ### Computing eta
    eta=zeta*recip_efr2

    ### Computing pressure pressure
    press = 0.5*(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+kt))*(1.-kt)

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def ellip_test(data={}, rheo={}, rheo_n = '', rot=False):
    '''
    ELLIPTICAL YIELD CURVE
    '''
    print('Computing Ellipse rheology test ; name:',rheo_n)


    # load data
    ep = data['ep']
    em = data['em']
    eI = data['eI']
    eII = data['eII']
    if rot :
        e12 = data['e12']
        e21 = data['e21']
    else:
        e12 = data['e12s']
        e21 = data['e12s']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'efr' in rheo:
        efr = rheo['efr']
    elif e != 2:
        efr = e
        rheo['efr'] = efr
    else:
        efr = e_d
        rheo['efr'] = efr

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    recip_e2 = 1./(e*e)
    recip_efr2 = 1./(efr*efr)
    e2_recip_efr4 = e**2/(efr**4)

    ### Computing Delta
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*np.abs(e12*e21))
    # deltaCsq=ep*ep+e2_recip_efr4*(np.abs(em*em+4.0*e12*e21))
    deltaCsq = ep**2 + recip_efr2 * (em**2 + 4.0*e12**2)
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*e12*e21)

    deltaC = np.sqrt(deltaCsq)
    ## with a sqrt function
    if SEAICE_DELTA_SMOOTHREG :
        deltaCreg = np.sqrt(deltaCsq + deltaMinSq)
    ## with a sharp min function
    else:
        deltaCreg = np.sqrt(np.maximum(deltaCsq,deltaMinSq))

    ### Computing Zeta
    ## with a smooth tanh function
    if SEAICE_ZETA_SMOOTHREG:
        argTmp = np.exp(-1./(deltaCreg*SEAICE_zetaMaxFac))
        zeta=ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+kt)
    ## with a sharp min function
    else:
        zeta=press0*(1.+kt)/(2.*deltaCreg)
        zeta=np.minimum(zeta,ZMAX)
    zeta=np.maximum(zeta,ZMIN)

    ### Computing pressure pressure
    press = 0.5*(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+kt))*(1.-kt)

    ### Computing eta
    sI = zeta * ep - press
    eta = 1. / e * np.sqrt( kt*press0**2 -  sI * ( sI + press0*(1 - kt) ) ) / (2 * eII + 1e-20)


    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def mce(data={}, rheo={}, rheo_n = ''):
    '''
    MC-E rheology
    '''
    print('Computing MC-E rheology ; name:',rheo_n)


    # load data
    eI = data['eI']
    eII  = data['eII']
    em = data['em']
    e12 = data['e12s']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    if 'mu_c' in rheo:
        mu_c = rheo['mu_c']
    else:
        mu_c = mu_c_d
        rheo['mu_c'] = mu_c

    recip_e2 = 1/e**2

    # compute delta
    deltaCsq = eI*eI+recip_e2*(np.abs(em*em+4.0*e12**2))

    deltaC = np.sqrt(deltaCsq)
    ## with a sqrt function
    if SEAICE_DELTA_SMOOTHREG :
        deltaCreg = np.sqrt(deltaCsq + deltaMinSq)
    ## with a sharp min function
    else:
        deltaCreg = np.sqrt(np.maximum(deltaCsq,deltaMinSq))

    ### Computing Zeta
    ## with a smooth tanh function
    if SEAICE_ZETA_SMOOTHREG:
        argTmp = np.exp(-1./(deltaCreg*SEAICE_zetaMaxFac))
        zeta = ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+kt)
    ## with a sharp min function
    else:
        zeta = press0*(1.+kt)/(2.*deltaCreg)
        zeta = np.minimum(zeta,ZMAX)
    zeta = np.maximum(zeta,ZMIN)

    ### Computing pressure pressure
    press = 0.5*(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+kt))*(1.-kt)

    ### Computing eta
    eta_mc = mu*(press-zeta*eI+kt*press0)/(np.maximum(deltaMin,eII))

    ### compressive capping eta
    eta_c = 0.5*mu_c*press0*(1+kt)*(eI/deltaCreg + 1.)/(np.maximum(deltaMin,eII))
    eta = np.minimum(eta_mc, eta_c)

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def mcs(data={}, rheo={}, rheo_n = ''):
    '''
    MC-S rheology
    '''
    print('Computing MC-S rheology ; name:',rheo_n)

    # load data
    eI = data['eI']
    em = data['em']
    e12 = data['e12s']
    eII = data['eII']

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    zetaMax = press0*(1+kt)/(2.*deltaMin)

    # compute the viscosities
    zeta = np.minimum(zetaMax,press0*(1+kt)/(np.fabs(2*eI+1e-20)))

    # press = (press0 * (1.-SEAICEpressReplFac) + SEAICEpressReplFac * 2 * zeta * np.fabs(ep)/(1.+kt))*(1.-kt)
    press = 0.5 * press0 * ( 1. - kt )

    eta = np.minimum(mu*(press-zeta*eI+kt*press0)/(np.maximum(1e-20,eII)),zetaMax)

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def mceG(data={}, rheo={}, rheo_n = ''):
    '''
    MC-EG rheology
    '''
    print('Computing MC-E rheology defined with a G function ; name:',rheo_n)

    # load data
    eI = data['eI']
    eII = data['eII']

    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    if 'mu_c' in rheo:
        mu_c = rheo['mu_c']
    else:
        mu_c = mu_c_d
        rheo['mu_c'] = mu_c

    # compute the MC zeta
    delta_mc = np.maximum(eII + mu * e**2 * eI, 1e-20)
    zeta_mc = (1. + kt) * mu * e**2 * press0 / (2. * delta_mc)

    # compute the cap zeta
    delta_c = np.maximum(eII - mu_c * e**2 * eI ,1e-20)
    zeta_c  = (1. + kt) * mu_c * e**2 * press0 / (2. * delta_c)

    # Apply viscous capping
    zeta_mc = np.minimum(zeta_mc,ZMAX)
    zeta_c  = np.minimum(zeta_c ,ZMAX)

    press_mc = 0.5 * (1. - kt) * ( press0 * (1. - SEAICEpressReplFac) + SEAICEpressReplFac * zeta_mc * 2 * delta_mc / ( (1+kt) * mu   * e**2))
    press_c =  0.5 * (1. - kt) * ( press0 * (1. - SEAICEpressReplFac) + SEAICEpressReplFac * zeta_c  * 2 * delta_c  / ( (1+kt) * mu_c * e**2))

    # compute MC and cap eta
    eta_c  = zeta_c  / e**2
    eta_mc = zeta_mc / e**2

    # Choose between MC limbs and cap
    eta = np.minimum(eta_mc, eta_c)
    zeta = np.minimum(zeta_mc, zeta_c)
    press = np.minimum(press_mc, press_c)

    # Choose between MC limbs and cap
    # eta = eta_mc #np.minimum(eta_mc, eta_c)
    # zeta = zeta_mc #np.minimum(zeta_mc, zeta_c)
    # press = press_mc #np.minimum(press_mc, press_c)

    # Choose between MC limbs and cap
    # eta = eta_c #np.minimum(eta_mc, eta_c)
    # zeta = zeta_c #np.minimum(zeta_mc, zeta_c)
    # press = press_c #np.minimum(press_mc, press_c)

    # Choose between MC limbs and cap
    # zeta = np.where(eta_mc<eta_c, zeta_mc, zeta_c)
    # press = np.where(eta_mc<eta_c, press_mc, press_c)
    # eta = np.where(eta_mc<eta_c, eta_mc, eta_c)

    ## When no cap is used, the Conditions to be in the triangle
    # lim_zeta = press0 / ( 4 * np.abs(eI) + 1e-20)
    # zeta = np.minimum(zeta,lim_zeta)
    # lim_eta = mu * press0 / ( 2*eII + 1e-20)
    # eta = np.minimum(eta,lim_eta)

    # press = zeta * delta_c * (1. - kt) / ((1.+kt) * e**2 * mu_c)
    # press = 0.5 * press0 * ( 1. - kt )

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None



def mctd(data={}, rheo={}, rheo_n = ''):
    '''
    MC-TD rheology
    '''
    print('Computing MC-TD rheology ; name:',rheo_n)


    # load data
    eI = data['eI']
    eII = data['eII']

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'mu_c' in rheo:
        mu_c = rheo['mu_c']
    else:
        mu_c = mu_c_d
        rheo['mu_c'] = mu_c

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    k = eI / ( eII + 1e-20)
    x = (-(6.*(1.+kt)-2*k*k)+2.*k*np.sqrt(k*k+3. * (1.+kt)))/9. + kt

    alpha = 0.99
    x = np.minimum( x, alpha*kt )
    # p = 100.
    # x = -np.log(np.exp(-p*x)+np.exp(-p*kt))/p

    cyc = (2. - kt)/3.

    zeta = ( x + cyc ) / np.copysign(np.maximum(np.fabs(eI), 1e-20),eI) * press0
    # zeta = 2 * press0 / (9 * eII + 1e-20) * (k + np.sqrt(k**2 + 3 * (1 - kt) ) )

    press = cyc * press0

    eta_mc = - press0 * mu * ( x - kt ) / ( eII + 1e-20 )
    # eta_mc = mu * (press - zeta * eI - press0 * kt) / ( 2 * eII + 1e-20 )

    # compressive cap
    eta_c = mu_c * (zeta * eI - press + press0) / (np.maximum(1e-20,eII))
    eta = np.minimum(eta_mc, eta_c)

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def mcpl(data={}, rheo={}, rheo_n = ''):
    '''
    MC-PL rheology
    '''
    print('Computing MC-PL rheology ; name:',rheo_n)


    # load data
    eI = data['eI']
    eII = data['eII']

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'mu_c' in rheo:
        mu_c = rheo['mu_c']
    else:
        mu_c = mu_c_d
        rheo['mu_c'] = mu_c


    k = eI / (eII + 1e-20) # 0

    x = 0.5 * (k - 1 + kt) #

    alpha=0.95
    x = np.minimum( x, alpha*kt )
    x = np.maximum( x, -1+(1-alpha)*kt )

    cyc = 0.5 * (1 - kt)

    zeta = ( x + cyc ) / np.copysign(np.maximum(np.fabs(eI), deltaMinSq),eI) * press0

    eta_mc = - press0 * mu * ( x - kt ) / ( eII + 1e-20 )

    press = cyc * press0

    eta_c = mu_c * (zeta * eI - press + press0) / (np.maximum(deltaMin,eII))

    eta = np.minimum(eta_mc, eta_c)

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def td(data={}, rheo={}, rheo_n = ''):
    '''
    TD rheology
    '''
    print('Computing TD rheology ; name:',rheo_n)

    # load data
    eI = data['eI']
    eII = data['eII']

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    k = eI / ( eII + 1e-20)
    x = (-(6.*(1.+kt)-2*k*k)+2.*k*np.sqrt(k*k+3. * (1.+kt)))/9. + kt

    alpha = 0.95

    x = np.minimum( x, alpha*kt )

    cyc = (2. - kt)/3.

    zeta = ( x + cyc ) / np.copysign(np.maximum(np.fabs(eI), deltaMinSq),eI) * press0

    eta = - ( x - kt ) * np.sqrt( 1. + x ) * press0 / (eII + 1e-20)

    press = cyc * press0

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def pl(data={}, rheo={}, rheo_n = ''):
    '''
    PL rheology
    '''
    print('Computing PL rheology ; name:',rheo_n)

    # load data
    eI = data['eI']
    eII = data['eII']

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0


    k = eI / (eII + 1e-20) # 0

    x = 0.5 * (k - 1 + kt) #

    alpha=0.95
    x = np.minimum( x, alpha*kt )
    x = np.maximum( x, -1+(1-alpha)*kt )

    cyc = 0.5 * (1 - kt)

    zeta = ( x + cyc ) / np.copysign(np.maximum(np.fabs(eI), 1e-20),eI) * press0

    eta = - ( x - kt ) * ( 1. + x ) * press0 / (eII + 1e-20)

    press = cyc * press0

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def etd(data={}, rheo={}, rheo_n = ''):
    '''
    E-TD rheology
    '''
    print('Computing E-TD rheology ; name:',rheo_n)

    # load data
    eI = data['eI']
    eII = data['eII']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    k = eI / ( eII + 1e-20)

    x = (-(6.*(1.+kt)-2*k*k)+2.*k*np.sqrt(k*k+3. * (1.+kt)))/9. + kt

    alpha = 0.99

    x = np.minimum( x, alpha*kt )

    cyc = (2. - kt) / 3.

    zeta = (x + cyc) / np.copysign(np.maximum(np.fabs(eI), 1e-20),eI) * press0

    eta = 1. / e * np.sqrt( kt - x * ( x + 1-kt ) ) / (eII + 1e-20) * press0

    press = cyc * press0

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def epl(data={}, rheo={}, rheo_n = ''):
    '''
    E-PL rheology
    '''
    print('Computing E-PL rheology ; name:',rheo_n)

    # load data
    eI = data['eI']
    eII = data['eII']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    k = eI / ( eII + 1e-20)

    x = 0.5 * (k - 1 + kt)

    alpha=0.99

    x = np.minimum( x, alpha*kt )
    x = np.maximum( x, -1+(1-alpha)*kt )

    cyc = 0.5 * (1 - kt)

    zeta = (x + cyc) / np.copysign(np.maximum(np.fabs(eI), 1e-20),eI) * press0

    eta=  1 / e * np.sqrt( kt - x * ( x + 1-kt ) ) / (eII + 1e-20) * press0

    press = cyc * press0

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def muID(data={}, rheo={}, rheo_n = ''):
    '''
    mu(I) rheological framework
    Heyman, J., Delannay, R., Tabuteau, H., & Valance, A. (2017). Compressibility regularizes the mu(I)-rheology for dense granular flows. Journal of Fluid Mechanics, 830, 553–568. https://doi.org/10.1017/jfm.2017.612
    Herman, A. (2022). Granular effects in sea ice rheology in the marginal ice zone. Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences, 380(2235), 20210260. https://doi.org/10.1098/rsta.2021.0260
    '''
    print('Computing mu(I) rheology ; name:',rheo_n)

    # load data
    eI = data['eI']
    eII = data['eII']
    phi = data['phi']
    h = data['h']

    # load rheo parameters
    if 'phi_0' in rheo:
        phi_0 = rheo['phi_0']
    else:
        phi_0 = phi_0_d
        rheo['phi_0'] = phi_0

    if 'c_phi' in rheo:
        c_phi = rheo['c_phi']
    else:
        c_phi = c_phi_d
        rheo['c_phi'] = c_phi

    if 'rho' in rheo:
        rho = rheo['rho']
    else:
        rho = rho_d
        rheo['rho'] = rho

    if 'd_m' in rheo:
        d_m = rheo['d_m']
    else:
        d_m = d_m_d
        rheo['d_m'] = d_m

    if 'mu_0' in rheo:
        mu_0 = rheo['mu_0']
    else:
        mu_0 = mu_0_d
        rheo['mu_0'] = mu_0

    if 'mu_i' in rheo:
        mu_i = rheo['mu_i']
    else:
        mu_i = mu_i_d
        rheo['mu_i'] = mu_i

    if 'mub_0' in rheo:
        mub_0 = rheo['mub_0']
    else:
        mub_0 = mub_0_d
        rheo['mub_0'] = mub_0

    if 'mub_i' in rheo:
        mub_i = rheo['mub_i']
    else:
        mub_i = mub_i_d
        rheo['mub_i'] = mub_i

    if 'I_0' in rheo:
        I_0 = rheo['I_0']
    else:
        I_0 = I_0_d
        rheo['I_0'] = I_0

    if 'coh' in rheo:
        coh = rheo['coh']
    else:
        coh = coh_d
        rheo['coh'] = coh

    if 'Pmax' in rheo:
        Pmax = rheo['Pmax']
    else:
        Pmax = Pmax_d
        rheo['Pmax'] = Pmax

    if 'Cstar' in rheo:
        Cstar = rheo['Cstar']
    else:
        Cstar = Cstar_d
        rheo['Cstar'] = Cstar

    # To force and constant phi
    # phi_d = 1.0

    # to print debug values
    if 'db' in rheo:
        db = rheo['db']
    else:
        db = False

    press0 = Pmax * h * np.exp(-Cstar*(1-phi))
    if db: pr_var_stats(press0, 'press0')

    I = (phi_0 - phi)/c_phi
    if db: pr_var_stats(I, 'I')

    # H22
    # alpha, beta are FSD parameters
    # def p(d, alpha, beta):
    #     return d**(-alpha) * np.exp(-d/beta)
    # c1 = fc1(alpha, beta) # unknown?
    # c2 = fc2(alpha, beta) # unknown?
    # I = ( 1. / c2 * np.arctanh(1. / c1 * (phi_0 - phi) ) )**2
    # def d_m(alpha, beta):
    #     return 2 * beta * (3 - alpha)

    muI = mu_0 + (mu_i-mu_0)/(I_0/I+1)
    if db: pr_var_stats(muI, 'muI ')

    mubI = mub_0 + (mub_i-mub_0)/(I_0/I+1)
    # mubI = 0
    if db: pr_var_stats(mubI, 'mubI ')

    # press = np.minimum(rho*(d*eII/(phi-phi_0))**2*press0,Pmax)
    press = rho*(d_m*eII/(phi-phi_0))**2 * press0
    # press = Pmax
    if db: pr_var_stats(press, 'press')

    p_eff = press * ( 1 - mubI * eI/ (eII + 1e-20) )
    if db: pr_var_stats(p_eff, 'p_eff')

    eta = press / (2 * eII+ 1e-20) * muI
    if db: pr_var_stats(eta, 'eta')

    zeta = press / (2 * eII+ 1e-20) * (muI + 2 * mubI)
    if db: pr_var_stats(zeta, 'zeta')

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press
    data[rheo_n]['press0'] = press0

    return None

##########
# STRESSES
##########

def compute_stress(data={}):

    for rheo_n in data['rheos']:
        compu_sigma(data=data, rheo_n=rheo_n)
        comp_str_inva(data=data, rheo_n=rheo_n)

    return None

def compu_sigma(data={}, rheo_n='', plot=False):

    zeta = data[rheo_n]['zeta']
    eta = data[rheo_n]['eta']
    press = data[rheo_n]['press']
    ep = data['ep']
    em = data['em']
    e21 = data['e21']
    e12 = data['e12']

    data[rheo_n]['sig11'] = zeta*ep+eta*em-press
    data[rheo_n]['sig22'] = zeta*ep-eta*em-press
    data[rheo_n]['sig12'] = 2*e21*eta
    data[rheo_n]['sig21'] = 2*e12*eta

    if plot :
        plt.figure('sig11')
        plt.pcolormesh(data[rheo_n]['sig11'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig22')
        plt.pcolormesh(data[rheo_n]['sig22'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig12')
        plt.pcolormesh(data[rheo_n]['sig12'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig21')
        plt.pcolormesh(data[rheo_n]['sig21'])
        plt.axis('equal')
        plt.colorbar()

        plt.show()

    return None

def comp_princ_stress(data={}, rheo_n=''):

    sig11 = data[rheo_n]['sig11']
    sig12 = data[rheo_n]['sig12']
    sig22 = data[rheo_n]['sig22']
    sig21 = data[rheo_n]['sig21']
    press0 = data[rheo_n]['press0']

    sigp=sig11+sig22
    sigm=sig11-sig22
    # sigTmp=np.sqrt(sigm**2+4*np.abs(sig12*sig21))
    sigTmp=np.sqrt(np.abs(sigm**2+4*sig12*sig21))
    # sigTmp=np.sqrt(sigm**2+4*sig12*sig21)

    data[rheo_n]['sig1']=0.5*(sigp+sigTmp)
    data[rheo_n]['sig2']=0.5*(sigp-sigTmp)

    data[rheo_n]['sig1n']=0.5*(sigp+sigTmp)/press0
    data[rheo_n]['sig2n']=0.5*(sigp-sigTmp)/press0

    return None

def comp_str_inva(data={}, rheo_n='', plot=False):

    comp_princ_stress(data=data, rheo_n=rheo_n)

    sig1 = data[rheo_n]['sig1']
    sig2 = data[rheo_n]['sig2']

    data[rheo_n]['sigI'] = 0.5*(sig1+sig2)
    data[rheo_n]['sigII'] = 0.5*(sig1-sig2)

    sig1n = data[rheo_n]['sig1n']
    sig2n = data[rheo_n]['sig2n']

    data[rheo_n]['sigIn'] = 0.5*(sig1n+sig2n)
    data[rheo_n]['sigIIn'] = 0.5*(sig1n-sig2n)

    if plot :
        plt.figure('sigI')
        plt.pcolormesh(data[rheo_n]['sigI'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sigII')
        plt.pcolormesh(data[rheo_n]['sigII'])
        plt.axis('equal')
        plt.colorbar()

        plt.show()

    return None

#########
# TESTING
#########

def compute_tests(data={}):

    for rheo_n in data['rheos']:
        for v in ['eta','zeta']:
            if (data[rheo_n][v] <= 0 ).any():
                print("WARNING: one (or more) viscosity", v, "is negative in ",rheo_n)

    return None

##########
# PLOTTING
##########

def plot_stress(data={}):

    fig1=plt.figure('stress states')
    ax = fig1.gca()
    plt.grid()
    plt.axis('equal')
    ax.set_ylabel('Sigma II (normalized)')
    ax.set_xlabel('Sigma I (normalized)')

    for rheo_n in data['rheos']:
        plot_inv(data=data, rheo_n=rheo_n, ax=ax, arrows=True)

    ax.legend(markerscale=5)
    ax.set_xlim([-1.2,0.2])

    return None

# def plot_inv(sigI,sigII,eI,eII,opt=None,arrows=False,ax=None,carg=None):
def plot_inv(data={}, rheo_n='', opt=None, arrows=False, ax=None, carg=None):
    '''
    Plotting the yield curve in invariant coordinate
    '''

    sigI = data[rheo_n]['sigIn']
    sigII = data[rheo_n]['sigIIn']
    eI = data['eI']
    eII = data['eII']

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if 'plot_inv' in data[rheo_n]:
        if data[rheo_n]['plot_inv'] :
            fac = -1
        else:
            fac = 1
    else:
        fac = 1


    if carg != None :
        ax.scatter(sigI.ravel(), fac*sigII.ravel(), carg, label=rheo_n)
    else:
        p = ax.plot(sigI.ravel(), fac*sigII.ravel(), '.', ms=1, label=rheo_n)
        carg = p[0].get_color()

    if arrows :
        qpfac=10
        eu=np.hypot(eI[::qpfac,::qpfac],fac*eII[::qpfac,::qpfac])
        ax.quiver(sigI[::qpfac,::qpfac],fac*sigII[::qpfac,::qpfac],eI[::qpfac,::qpfac]/eu,fac*eII[::qpfac,::qpfac]/eu, scale=10, color=carg)

    if opt=='ellipse':
        t=tnsFac
        f=SEAICE_strength
        f=1
        elli=Ellipse(xy=((-f+t)/2.,0), width=(f+t)/e, height=(f+t), angle=-90,edgecolor='b', fc='None', lw=0.5)
        ax.add_patch(elli)

    if data[rheo_n]['rheo_t'] in ['mceG','ell']:
        if data[rheo_n]['rheo_t'] == 'mceG' : c = '-b'
        if data[rheo_n]['rheo_t'] == 'ell'  :
            if data[rheo_n]['e'] == data[rheo_n]['efr']:
                c = '-r'
            else:
                c = '-b'
        qpfac=10
        sI = np.array(sigI[::qpfac,::qpfac])
        sII = fac*np.array(sigII[::qpfac,::qpfac])
        kt = data[rheo_n]['kt']
        cx = -0.5 * (1 - kt) * np.ones(np.shape(sI))
        cy = 0 * np.ones(np.shape(sI))
        xs = np.stack((cx.flatten(),sI.flatten()),axis=1)
        ys = np.stack((cy.flatten(),sII.flatten()),axis=1)
        for i in range(len(xs)):
            ax.plot(xs[i], ys[i], c, lw = 1, alpha=0.3)

    return None


def plot_FR(data={}):

    fig1=plt.figure('flow rule')
    ax = fig1.gca()
    plt.grid()
    ax.set_ylabel('dilatancy angle delta=arctan(eI/eII) [$^\circ$]')
    ax.set_xlabel('Sigma I (normalized)')

    for rheo_n in data['rheos']:
        print('rheo_n FR ',rheo_n)
        plot_sIFR(data=data, rheo_n=rheo_n, ax=ax)

    ax.legend(markerscale=2)
    ax.set_ylim([-90,90])

    return None

def plot_sIFR(data={}, rheo_n='', ax=None, carg=None, opt=None):

    '''
    Plotting the the flow rule as function of sI
    '''

    sigI = data[rheo_n]['sigIn']
    eI = data['eI']
    eII = data['eII']

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if carg != None :
        ax.plot(sigI.ravel(),np.arctan((eI/eII).ravel())*180/np.pi,'.', color=carg, ms=4, label=rheo_n, alpha=0.2)
        # ax.plot(np.arctan((eI/eII).ravel())*180/np.pi,sigI.ravel(),'.', color=carg, ms=4, label=rheo_n, alpha=0.2)
    else:
        p = ax.plot(sigI.ravel(),np.arctan((eI/eII).ravel())*180/np.pi,'.', ms=4, label=rheo_n, alpha=0.2)
        # p = ax.plot(np.arctan((eI/eII).ravel())*180/np.pi,sigI.ravel(),'.', ms=4, label=rheo_n, alpha=0.2)
        carg = p[0].get_color()

    if data[rheo_n]['rheo_t'] == 'ell' :
        e = data[rheo_n]['e']
        efr = data[rheo_n]['efr']
        if efr != e :
            press = data[rheo_n]['press0']
            kt = data[rheo_n]['kt']
            fr = fr_th_ell(sigI, e, efr, press, kt)
            # ax.plot(np.arctan(fr.ravel())*180/np.pi, sigI.ravel(), 'xk', ms=6, label='th_nnfr_ell', alpha=0.3)
            ax.plot(sigI.ravel(),np.arctan(fr.ravel())*180/np.pi, 'xk', ms=2, label='th_nnfr_ell', alpha=0.2)

    return None


def fr_th_ell(sI, e, efr, press, kt):
    fr = e * (sI + 0.5 * (1 - kt)) / (efr**2 * np.sqrt(kt - sI * (sI + 1 - kt) ) )
    return fr

def plot_prAng(data={}):


    fig1=plt.figure('principal axes orientation')
    ax = fig1.gca()
    plt.grid()
    # plt.axis('equal')
    ax.set_xlabel('principal stress orientation')
    ax.set_ylabel('principal strain rate orientation')

    for rheo_n in data['rheos']:
        plot_prAng_ori(data=data, rheo_n=rheo_n, ax=ax)

    ax.legend(markerscale=2)
    ax.set_aspect('equal')

    return None

def plot_prAng_ori(data={}, rheo_n='', ax=None, carg=None, opt=None):

    '''
    Plotting the the flow rule as function of sI
    '''

    sig11 = data[rheo_n]['sig11']
    sig22 = data[rheo_n]['sig22']
    sig12 = data[rheo_n]['sig12']

    e11 = data['e11']
    e22 = data['e22']
    e12 = data['e12']

    psi_st = 0.5 * np.arctan2(2*sig12,(sig11-sig22)) * 180/np.pi
    psi_sr = 0.5 * np.arctan2(2*e12,(e11-e22)) * 180/np.pi

    # psi_st = 0.5 * np.arctan(2*sig12/(sig11-sig22)) * 180/np.pi
    # psi_sr = 0.5 * np.arctan(2*e12/(e11-e22)) * 180/np.pi

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if carg != None :
        ax.plot(psi_st.ravel(),psi_sr.ravel()-psi_st.ravel(),'.', color=carg, ms=4, label=rheo_n, alpha=0.2)
        # ax.plot(psi_st.ravel(),psi_sr.ravel(),'.', color=carg, ms=4, label=rheo_n, alpha=0.2)
    else:
        # p = ax.plot(psi_st.ravel(),psi_sr.ravel()-psi_st.ravel(),'.', ms=4, label=rheo_n, alpha=0.2)
        p = ax.plot(psi_st.ravel(),psi_sr.ravel(),'.', ms=4, label=rheo_n, alpha=0.2)
        carg = p[0].get_color()

    return None

####################
# OTHER RANDOM TOOLS
####################

def TD(sI,p=1):
    '''
    Shape of the TD yield curve
    - sII/P as function as sI/P
    '''
    return -(sI/p-a)*(1+sI/p)**q

def pr_var_stats(var, varname):
    print(varname, ' mean: ', np.nanmean(var), ' std: ', np.nanstd(var), ' min: ', np.nanmin(var), ' max: ', np.nanmax(var))


