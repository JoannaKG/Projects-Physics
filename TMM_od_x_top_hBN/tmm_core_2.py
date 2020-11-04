from __future__ import division, print_function, absolute_import

from numpy import cos, inf, zeros, array, exp, conj, nan, isnan, pi, sin

import numpy as np
import scipy as sp
from scipy.integrate import quad
import csv
from scipy.interpolate import interp1d
degree = pi/180

import sys
EPSILON = sys.float_info.epsilon # typical floating-point calculation error

def make_2x2_array(a, b, c, d, dtype=float):
    """
    Makes a 2x2 numpy array of [[a,b],[c,d]]

    Same as "numpy.array([[a,b],[c,d]], dtype=float)", but ten times faster
    """
    my_array = np.empty((2,2), dtype=dtype)
    my_array[0,0] = a
    my_array[0,1] = b
    my_array[1,0] = c
    my_array[1,1] = d
    return my_array

def is_forward_angle(n, theta):
    """
    if a wave is traveling at angle theta from normal in a medium with index n,
    calculate whether or not this is the forward-traveling wave (i.e., the one
    going from front to back of the stack, like the incoming or outgoing waves,
    but unlike the reflected wave). For real n & theta, the criterion is simply
    -pi/2 < theta < pi/2, but for complex n & theta, it's more complicated.
    See https://arxiv.org/abs/1603.02720 appendix D. If theta is the forward
    angle, then (pi-theta) is the backward angle and vice-versa.
    """
    assert n.real * n.imag >= 0, ("For materials with gain, it's ambiguous which "
                                  "beam is incoming vs outgoing. See "
                                  "https://arxiv.org/abs/1603.02720 Appendix C.\n"
                                  "n: " + str(n) + "   angle: " + str(theta))
    ncostheta = n * cos(theta)
    if abs(ncostheta.imag) > 100 * EPSILON:
        # Either evanescent decay or lossy medium. Either way, the one that
        # decays is the forward-moving wave
        answer = (ncostheta.imag > 0)
    else:
        answer = (ncostheta.real > 0)
    answer = bool(answer)
    error_string = ("It's not clear which beam is incoming vs outgoing. Weird"
                    " index maybe?\n"
                    "n: " + str(n) + "   angle: " + str(theta))
    if answer is True:
        assert ncostheta.imag > -100 * EPSILON, error_string
        assert ncostheta.real > -100 * EPSILON, error_string
        assert (n * cos(theta.conjugate())).real > -100 * EPSILON, error_string
    else:
        assert ncostheta.imag < 100 * EPSILON, error_string
        assert ncostheta.real < 100 * EPSILON, error_string
        assert (n * cos(theta.conjugate())).real < 100 * EPSILON, error_string
    return answer

def snell(n_1, n_2, th_1):
    # Important that the arcsin here is scipy.arcsin, not numpy.arcsin! (They
    # give different results e.g. for arcsin(2).)
    th_2_guess = sp.arcsin(n_1*np.sin(th_1) / n_2)
    if is_forward_angle(n_2, th_2_guess):
        return th_2_guess
    else:
        return pi - th_2_guess

def list_snell(n_list, th_0):
    angles = sp.arcsin(n_list[0]*np.sin(th_0) / n_list)
    # The first and last entry need to be the forward angle (the intermediate
    # layers don't matter, see https://arxiv.org/abs/1603.02720 Section 5)
    if not is_forward_angle(n_list[0], angles[0]):
        angles[0] = pi - angles[0]
    if not is_forward_angle(n_list[-1], angles[-1]):
        angles[-1] = pi - angles[-1]
    return angles

def interface_r(polarization, n_i, n_f, th_i, th_f):

    if polarization == 's':
        return ((n_i * cos(th_i) - n_f * cos(th_f)) /
                (n_i * cos(th_i) + n_f * cos(th_f)))
    elif polarization == 'p':
        return ((n_f * cos(th_i) - n_i * cos(th_f)) /
                (n_f * cos(th_i) + n_i * cos(th_f)))
    else:
        raise ValueError("Polarization must be 's' or 'p'")

def interface_t(polarization, n_i, n_f, th_i, th_f):

    if polarization == 's':
        return 2 * n_i * cos(th_i) / (n_i * cos(th_i) + n_f * cos(th_f))
    elif polarization == 'p':
        return 2 * n_i * cos(th_i) / (n_f * cos(th_i) + n_i * cos(th_f))
    else:
        raise ValueError("Polarization must be 's' or 'p'")

def R_from_r(r):

    return abs(r)**2

def T_from_t(pol, t, n_i, n_f, th_i, th_f):

    if pol == 's':
        return abs(t**2) * (((n_f*cos(th_f)).real) / (n_i*cos(th_i)).real)
    elif pol == 'p':
        return abs(t**2) * (((n_f*conj(cos(th_f))).real) /
                                (n_i*conj(cos(th_i))).real)
    else:
        raise ValueError("Polarization must be 's' or 'p'")


def interface_R(polarization, n_i, n_f, th_i, th_f):
    """
    Fraction of light intensity reflected at an interface.
    """
    r = interface_r(polarization, n_i, n_f, th_i, th_f)
    return R_from_r(r)

def interface_T(polarization, n_i, n_f, th_i, th_f):
    """
    Fraction of light intensity transmitted at an interface.
    """
    t = interface_t(polarization, n_i, n_f, th_i, th_f)
    return T_from_t(polarization, t, n_i, n_f, th_i, th_f)

def import_n(layer_csv):
    wavelength = []
    n = []
    k = []
    
    with open(layer_csv) as op:
        reader = csv.reader(op, delimiter=";")
        wavelength_list = list(zip(*reader))[0]
    for line in wavelength_list:
        wavelength.append(float(line)*1000)
        
    with open(layer_csv) as op:
        reader = csv.reader(op, delimiter=";")
        n_list = list(zip(*reader))[1]
    for line in n_list:
        n.append(float(line))
        
    with open(layer_csv) as op:
        reader = csv.reader(op, delimiter=";")
        k_list = list(zip(*reader))[2]
    for line in k_list:
        k.append(float(line)*1j)
    #print(wavelength, n, k)
        
    wavelength_array = array(wavelength)
    n_array = array(n)
    k_array = array(k)
        
    layer_n_list = []
    for i in range (0, len(wavelength)-1):
        layer_n_list.append([wavelength_array[i],n_array[i]+k_array[i]])

    layer_n_data = array(layer_n_list)
    layer_n_fn = interp1d(layer_n_data[:,0], layer_n_data[:,1], kind='linear')
    return layer_n_fn

def coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0=514):
    
    n_list = array(n_list)
    n_list_abs=array(n_list_abs)
    d_list = array(d_list, dtype=float)
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        n_layer = import_n('WSe2.csv')

    # Input tests
    if ((hasattr(lam_vac, 'size') and lam_vac.size > 1)
          or (hasattr(th_0, 'size') and th_0.size > 1)):
        raise ValueError('This function is not vectorized; you need to run one '
                         'calculation at a time (1 wavelength, 1 angle, etc.)')
    if (n_list.ndim != 1) or (d_list.ndim != 1) or (n_list.size != d_list.size):
        raise ValueError("Problem with n_list or d_list!")
    assert d_list[0] == d_list[-1] == inf, 'd_list must start and end with inf!'
    assert abs((n_list[0]*np.sin(th_0)).imag) < 100*EPSILON, 'Error in n0 or th0!'
    assert is_forward_angle(n_list[0], th_0), 'Error in n0 or th0!'
    num_layers = n_list.size

    # th_list is a list with, for each layer, the angle that the light travels
    # through the layer. Computed with Snell's law. Note that the "angles" may be
    # complex!
    th_list = list_snell(n_list, th_0)
    th_list_abs = list_snell(n_list_abs, th_0)

    # kz is the z-component of (complex) angular wavevector for forward-moving
    # wave. Positive imaginary part means decaying.
    kz_list = 2 * np.pi * n_list * cos(th_list) / lam_vac
    kz_list_abs = 2 * np.pi * n_list_abs * cos(th_list_abs) / lam0


    olderr = sp.seterr(invalid='ignore')
    delta = kz_list * d_list
    delta_abs = kz_list_abs * d_list
    
    def beta_x(x, lam, dtype=complex):
        beta_x=-2*np.pi*TMD_n_fn(lam)*x/lam
        return beta_x
    
    sp.seterr(**olderr)

    for i in range(1, num_layers-1):
        if delta[i].imag > 35:
            delta[i] = delta[i].real + 35j
            if 'opacity_warning' not in globals():
                global opacity_warning
                opacity_warning = True
                print("Warning: Layers that are almost perfectly opaque "
                      "are modified to be slightly transmissive, "
                      "allowing 1 photon in 10^30 to pass through. It's "
                      "for numerical stability. This warning will not "
                      "be shown again.")
                
    #--------------liczenie Fsc(lam_vac)-------------------------------            

    t_list = zeros((num_layers, num_layers), dtype=complex)
    r_list = zeros((num_layers, num_layers), dtype=complex)
    for i in range(num_layers-1):
        t_list[i,i+1] = interface_t(pol, n_list[i], n_list[i+1],
                                    th_list[i], th_list[i+1])
        #print("t_list["+str(i)+","+str(i+1)+"] ="+str(t_list[i,i+1]))
        r_list[i,i+1] = interface_r(pol, n_list[i], n_list[i+1],
                                    th_list[i], th_list[i+1])
        #print("r_list["+str(i)+","+str(i+1)+"] ="+str(r_list[i,i+1]))
        
    t_list[1,0]=interface_t(pol, n_list[1], n_list[0],th_list[1], th_list[0])
        
    M_list = zeros((num_layers, 2, 2), dtype=complex)
    for i in range(1, num_layers-1):
        M_list[i] = (1/t_list[i,i+1]) * np.dot(
            make_2x2_array(exp(-1j*delta[i]), 0, 0, exp(1j*delta[i]),
                           dtype=complex),
            make_2x2_array(1, r_list[i,i+1], r_list[i,i+1], 1, dtype=complex))
    Mtilde = make_2x2_array(1, 0, 0, 1, dtype=complex)
    Mtilde1 = make_2x2_array(1, 0, 0, 1, dtype=complex)
    
    for i in range(1, num_layers-1):
        Mtilde = np.dot(Mtilde, M_list[i])
        Mtilde1 = np.dot(Mtilde1, M_list[i])
    Mtilde = np.dot(make_2x2_array(1, r_list[0,1], r_list[0,1], 1,
                                   dtype=complex), Mtilde)
   
    def F_sc(x,Mtilde1, dtype=complex):
        F_sc_arr = zeros((num_layers, 2, 2), dtype=complex)
        #print("t[1,0]="+str(t_list[1,0]))
        F_sc_arr = np.dot(make_2x2_array(1, r_list[0,1], t_list[1,0]*exp(1j*beta_x(x,lam_vac)), t_list[1,0]*exp(-1j*beta_x(x, lam_vac)),
                                   dtype=complex), Mtilde1)
        F_sc=F_sc_arr[1,0]/F_sc_arr[0,0]
        return F_sc
    
    #-------------liczenie Fabs(lam_0)------------------------------------
    t_list_abs = zeros((num_layers, num_layers), dtype=complex)
    r_list_abs = zeros((num_layers, num_layers), dtype=complex)
    for i in range(num_layers-1):
        t_list_abs[i,i+1] = interface_t(pol, n_list_abs[i], n_list_abs[i+1],
                                    th_list_abs[i], th_list_abs[i+1])
        #print("t_list["+str(i)+","+str(i+1)+"] ="+str(t_list[i,i+1]))
        r_list_abs[i,i+1] = interface_r(pol, n_list_abs[i], n_list_abs[i+1],
                                    th_list_abs[i], th_list_abs[i+1])
        #print("r_list["+str(i)+","+str(i+1)+"] ="+str(r_list[i,i+1]))
        
    t_list_abs[1,0]=interface_t(pol, n_list_abs[1], n_list_abs[0],th_list_abs[1], th_list_abs[0])
        
    M_list_abs = zeros((num_layers, 2, 2), dtype=complex)
    for i in range(1, num_layers-1):
        M_list_abs[i] = (1/t_list_abs[i,i+1]) * np.dot(
            make_2x2_array(exp(-1j*delta_abs[i]), 0, 0, exp(1j*delta_abs[i]),
                           dtype=complex),
            make_2x2_array(1, r_list_abs[i,i+1], r_list_abs[i,i+1], 1, dtype=complex))
    Mtilde_abs = make_2x2_array(1, 0, 0, 1, dtype=complex)
    
    for i in range(1, num_layers-1):
        Mtilde_abs = np.dot(Mtilde_abs, M_list_abs[i])

    def F_abs(x,Mtilde_abs, dtype=complex):
        F_abs_arr = zeros((num_layers, 2, 2), dtype=complex)
        F_abs_arr = np.dot(make_2x2_array(1, r_list_abs[0,1], t_list_abs[0,1]*exp(1j*beta_x(x, lam0)), t_list_abs[0,1]*exp(-1j*beta_x(x, lam0)),
                                   dtype=complex), Mtilde_abs)
        F_abs=F_abs_arr[1,0]/F_abs_arr[0,0]
        
        return F_abs
    
    
    def integrand(x):
            return (abs(F_abs(x,Mtilde_abs)*F_sc(x,Mtilde1)))**2    
            #return (abs(F_abs(x,Mtilde_abs)))**2    
    ans, err = quad(integrand,0, layer_thickness)
    Enh=ans
    #print("lam_vac="+str(lam_vac)+", integrand="+str(ans))

    r = Mtilde[1,0]/Mtilde[0,0]
    t = 1/Mtilde[0,0]

    R = R_from_r(r)

    T = T_from_t(pol, t, n_list[0], n_list[-1], th_0, th_list[-1])
 

    return {'r': r, 't': t, 'R': R, 'T': T, 'Enh' : Enh, 'kz_list': kz_list, 'th_list': th_list,
            'pol': pol, 'n_list': n_list, 'd_list': d_list, 'th_0': th_0,
            'lam_vac':lam_vac}