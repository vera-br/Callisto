# induced field

import numpy as np
from scipy import constants
import scipy.special as sps
from maths import *
from field_functions import *
from Aeiphi import *

mu0 = constants.mu_0
pi = constants.pi

def R_bessel(r1, k):
    z = r1 * k
    J_min_52 = sps.jv(-5 / 2, z)
    J_32 = sps.jv(3 / 2, z)
    J_12 = sps.jv(1 / 2, z)

    return (z * J_min_52) / (3 * J_32 - z * J_12)


def ae_iphi(sigma, omega, rm, r0, r1):
    k = (1 - 1j) * np.sqrt(mu0 * sigma * omega / 2)
    R = R_bessel(r1, k)

    amp = (r0 / rm) ** 3
    z = r0 * k

    J_52 = sps.jv(5 / 2, z)
    J_min_52 = sps.jv(-5 / 2, z)
    J_12 = sps.jv(1 / 2, z)
    J_min_12 = sps.jv(-1 / 2, z)

    return np.real(amp * (R * J_52 - J_min_52) / (R * J_12 - J_min_12))


def B_induced_finite_conductivity(orbit, B_external, sigma, omega, Rm, R0, R1):
    """
    Calculate the induced magnetic field with finite conductivity
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param sigma: conductivity in S
    :param omega: angular frequency of inducing field in
    :param rm: object radius in m
    :param r0: conducting layer outer radius in m
    :param r1: conducting layer inner radius in m
    :return: time evolution array of Bx, By, Bz in nT
    for more info see (Zimmer et al) https://www.sciencedirect.com/science/article/pii/S001910350096456X
    """
    orbit = orbit.transpose()

    A = ae_iphi(sigma, omega, Rm, R0, R1)

    Bind_evolution = []
    for B_ext, vector in zip(B_external, orbit):

        position = vector[1:4] 

        M = -(2 * pi / mu0) * A * B_ext * (Rm**3)

        rmag = np.linalg.norm(position)
        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)

def B_induced_finite_conductivity2(orbit, B_external, ocean_sigma, omega, ocean_depth, surface_layer):
    """
    Calculate the induced magnetic field with finite conductivity
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param sigma: conductivity in S
    :param omega: angular frequency of inducing field in
    :param rm: object radius in m
    :param r0: conducting layer outer radius in m
    :param r1: conducting layer inner radius in m
    :return: time evolution array of Bx, By, Bz in nT
    for more info see (Zimmer et al) https://www.sciencedirect.com/science/article/pii/S001910350096456X
    """
    orbit = orbit.transpose()

    A = Aeiphi(ocean_depth, surface_layer, ocean_sigma, omega)

    Bind_evolution = []
    for B_ext, vector in zip(B_external, orbit):

        position = vector[1:4] 

        M = -(2 * pi / mu0) * A * B_ext * (R_C**3)

        rmag = np.linalg.norm(position)
        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)


def B_induced_infinite(orbit, B_external, Rm, R0):
    """
    Induced magnetic field for infinite conductivity (superconductor)
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param B_external: array of external field vectors Bx, By, Bz in nT
    :param Rm: object radius in m
    :param R0: conducting layer outer radius in m
    :return: time evolution array of Bx, By, Bz in nT
    """
    orbit = orbit.transpose()

    A = -((R0 / Rm) ** 3)

    B_ind_evolution = []
    for B_ext, vector in zip(B_external, orbit):

        pos = vector[1:4] 

        M = -(2 * pi / mu0) * A * B_ext * (Rm**3)

        rmag = np.linalg.norm(pos)
        rdotM_r = np.dot(pos, M) * pos

        B_ind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        B_ind_evolution.append(B_ind)

    return np.array(B_ind_evolution)

def ae_iphi_multilayer(conductivities, r, l, omega):
    '''
    Calculates Ae^iphi for a given l and omega based on the approach given by Seufert et al. (2011)
    '''
    # def Fs(l, r, k):
    #     rk = r * k
    #     sqrt = np.sqrt(pi / (2 * rk))
    #     I_v = sps.iv(l + 0.5, rk)
    #     dI_v = sps.ivp(l + 0.5, rk)
    #     K_v = sps.kv(l + 0.5, rk)
    #     dK_v = sps.kvp(l + 0.5, rk)
    #     F1 = sqrt * I_v
    #     F2 = sqrt * K_v
    #     dF1 = -0.5 * F1 / r + sqrt * dI_v
    #     dF2 = -0.5 * F2 / r + sqrt * dK_v
    #     return F1, F2, dF1, dF2
    
    def Fs(l, r, k):
        rk = r * k
        F1 = sps.spherical_in(l, rk)
        F2 = sps.spherical_kn(l, rk)
        dF1 = sps.spherical_in(l, rk, derivative=True)
        dF2 = sps.spherical_in(l, rk, derivative=True)
        return F1, F2, dF1, dF2

    k1 = np.sqrt(-1j * omega * mu0 * conductivities[0])
    k2 = np.sqrt(-1j * omega * mu0 * conductivities[1])

    # Zimmer k
    # k1 = (1-1j) * np.sqrt(omega * mu0 * conductivities[0] / 2)
    # k2 = (1-1j) * np.sqrt(omega * mu0 * conductivities[1] / 2)
    
    r1 = r[0]
    F11, F21, dF11, dF21 = Fs(l, r1, k1)
    F12, F22, dF12, dF22 = Fs(l, r1, k2)

    Djmin1_Cjmin1 = (F12 / F22) * ( ( (dF12 / F12) - (dF11 / F11) ) / ( (dF11 / F11) - (dF22 / F22)) )
    
    for i in range(1, len(r) - 1):
       
       k_jmin1 = np.sqrt(-1j * omega * mu0 * conductivities[i])
       k_j = np.sqrt(-1j * omega * mu0 * conductivities[i + 1])

       # Zimmer k
    #    k_jmin1 = (1-1j) * np.sqrt(omega * mu0 * conductivities[i] / 2)
    #    k_j = (1-1j) * np.sqrt(omega * mu0 * conductivities[i + 1] / 2)

       r_jmin1 = r[i]

       F11, F21, dF11, dF21 = Fs(l, r_jmin1, k_jmin1)
       F12, F22, dF12, dF22 = Fs(l, r_jmin1, k_j)

    #    numer = (dF12 / F12) - (dF11 / F11) + Djmin1_Cjmin1 * (F21 / F11) * ( (dF12 / F12) - (dF21 / F21))
    #    denom = (dF11 / F11) - (dF22 / F22) + Djmin1_Cjmin1 * (F21 / F11) * ( (dF21 / F21) - (dF22 / F22))
    #    Djmin1_Cjmin1 = (F12 / F22) * (numer / denom)
       
       Djmin1_Cjmin1 = F12 / F22 * ((dF12 / F12) - (dF11 / F11) + Djmin1_Cjmin1 * (F21 / F11) * ( (dF12 / F12) - (dF21 / F21))) / ((dF11 / F11) - (dF22 / F22) + Djmin1_Cjmin1 * (F21 / F11) * ( (dF21 / F21) - (dF22 / F22)))
    
    R = r[-1]

    k_J = np.sqrt(-1j * omega * mu0 * conductivities[-1])
    # k_J = (1-1j) * np.sqrt(omega * mu0 * conductivities[-1] / 2)

    F1J, F2J, dF1J, dF2J = Fs(l, R, k_J)

    # numer = (dF1J / F1J) - (l + 1) + Djmin1_Cjmin1 * (F2J / F1J) * ( (dF2J / F2J) - (l + 1) )
    # denom = (dF1J / F1J) + l + Djmin1_Cjmin1 * (F2J / F1J) * ( (dF2J / F2J) + l )
    # Ae_iphi = numer / denom
    Ae_iphi = ((dF1J / F1J) - (l + 1) + Djmin1_Cjmin1 * (F2J / F1J) * ( (dF2J / F2J) - (l + 1) )) / ((dF1J / F1J) + l + Djmin1_Cjmin1 * (F2J / F1J) * ( (dF2J / F2J) + l ))
    
    # print('|A| = {}'.format(abs(Ae_iphi)))
    # print('phi = {}'.format(np.arctan(Ae_iphi.imag / Ae_iphi.real) * 180 / np.pi))
    return Ae_iphi

def B_induced_finite_conductivity_multilayer(O, B_external, omega, conductivities, radii, Styczinski=False, aeiphi=None, shifted=False):
    """
    Calculate the induced magnetic field with finite conductivity
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param omega: angular frequency of inducing field in
    :param conductivities: array of conductivities of the layers in S
    :param radii: array of radii of the layers in m
    :return: time evolution array of Bx, By, Bz in nT
    """
    orbit = O.copy()
    t = orbit[0]
    orbit = orbit.transpose()
    if aeiphi != None:
        A = aeiphi
    elif Styczinski !=False:
        if Styczinski == True:
            A = Aeiphi_Styczinski(conductivities, radii, 1, omega)
        elif Styczinski == 'many':
            A = Aeiphi_Styczinski_many(conductivities, radii, 1, omega)
    else:
        A = ae_iphi_multilayer(conductivities, radii, 1, omega)
    
    B_ext = B_external.copy()
    B_ext[:,2] = 0
    if shifted == True:
        phi = -np.angle(A)
        t_phi = phi / omega
        t_interval = t[1] - t[0]
        n_intervals = int(np.round(t_phi / t_interval))
        B_ext = np.roll(B_ext, n_intervals, axis=0)

    _M = M = -(2 * pi / mu0) * abs(A) * (R_C**3)
    
    # _M = M = -(2 * pi / mu0) * A * (radii[-1]**3)
    
    # print('M = {}'.format(_M))

    Bind_evolution = []
    for B_ext, vector in zip(B_ext, orbit):

        position = vector[1:4]
        # M = -(2 * pi / mu0) * A * B_ext * (radii[-1]**3)
        M = B_ext * _M

        rmag = np.linalg.norm(position)
        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind = Bind.real
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)

def B_induced_finite_conductivity_multilayer_G(O, B_external, omega, conductivities, radii, Styczinski=False, aeiphi=None, shifted=False, t_longperiod=None):
    """
    Calculate the induced magnetic field with finite conductivity
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param omega: angular frequency of inducing field in
    :param conductivities: array of conductivities of the layers in S
    :param radii: array of radii of the layers in m
    :return: time evolution array of Bx, By, Bz in nT
    """
    orbit = O.copy()
    B_ext = B_external.copy()
    if shifted == True:
        B_ext[:,2] = B_ext[:,2] - np.mean(B_ext[:,2].ravel())
    else:
        B_ext[:,2] = 0
    t_red = orbit[0]
    
    orbit = orbit.transpose()
    if aeiphi != None:
        A = aeiphi
    elif Styczinski !=False:
        if Styczinski == True:
            A = Aeiphi_Styczinski(conductivities, radii, 1, omega)
        elif Styczinski == 'many':
            A = Aeiphi_Styczinski_many(conductivities, radii, 1, omega)
    else:
        A = ae_iphi_multilayer(conductivities, radii, 1, omega)
    
    if np.any(t_longperiod) != None:
        t = t_longperiod
        if shifted == True:
            # print(A)
            phi = -np.angle(A)
            # print(phi)
            # print(omega)
            t_phi = phi / omega
            # print(t_phi)
            t_interval = t[1] - t[0]
            # print(t_interval)
            try:
                n_intervals = int(np.round(t_phi / t_interval))
                factor = 1
            except:
                n_intervals = 0
                factor = 0
            B_ext = np.roll(B_ext, n_intervals, axis=0)

        tB_ext = np.c_[t_longperiod, B_ext]
        tB_ext = tB_ext.transpose()

        tB_reduced = []
        for i in range(len(t_red)):
            index = find_nearest_index(tB_ext[0], t_red[i])
            tB_i = tB_ext[:, int(index)]
            tB_reduced.append(tB_i)
        tB_reduced = np.transpose(tB_reduced)
        B_reduced = tB_reduced[1:]
        B_ext = B_reduced.transpose()

    # _M = M = -(2 * pi / mu0) * A * (radii[-1]**3)
    _M = M = -(2 * pi / mu0) * abs(A) * (R_C**3)
    # print('M = {}'.format(_M))

    Bind_evolution = []
    for B_ext, vector in zip(B_ext, orbit):

        position = vector[1:4]
        # M = -(2 * pi / mu0) * A * B_ext * (radii[-1]**3)
        M = B_ext * _M

        rmag = np.linalg.norm(position)
        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind = Bind.real
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution) * factor

def Aeiphi_Styczinski(conductivities, rs, n, omega):
    def get_spherical_harmonics(kr, n, derivatives=False):
        jn = sps.spherical_jn(n, kr)
        yn = sps.spherical_yn(n, kr)
        if derivatives == True:
            d_jn = sps.spherical_jn(n, kr, derivative=True)
            d_yn = sps.spherical_yn(n, kr, derivative=True)
            return jn, d_jn, yn, d_yn
        return jn, yn

    # inner boundary
    k1 = np.sqrt(1j * omega * mu0 * conductivities[0])
    r1 = rs[0]
    k2 = np.sqrt(1j * omega * mu0 * conductivities[1])

    cutoff_factor = 10
    skin_depth_factor = 10
    
    if np.abs(k1 * r1) > n * cutoff_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n, derivatives=True)
        del_l = -(jul + djul / (1j * k1 * r1)) / (yul + dyul / (1j * k1 * r1))
        print('Layer 1: kr >> n')
        if 1 / k1.imag > r1 / skin_depth_factor and 1 / k1.imag < r1 * skin_depth_factor:
            print('Approx. Invalid - 1/Im(k) = {}, layer thickness = {}'.format(1/k1.imag, r1))

    elif np.abs(k1 * r1) < n / cutoff_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        del_l = -(jul / yul)
        print('Layer 1: kr << n')
        if 1 / k1.imag > r1 / skin_depth_factor and 1 / k1.imag < r1 * skin_depth_factor:
            print('Approx. Invalid - 1/Im(k) = {}, layer thickness = {}'.format(1/k1.imag, r1))

    else:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
        delta_ul = jul * djll - jll * djul
        beta_ul = jll * dyul - yul * djll
        del_l = delta_ul / beta_ul

    for i in range(1, len(rs) - 1):
        kj = np.sqrt(1j * omega * mu0 * conductivities[i])
        rj = rs[i]

        rl = rs[i-1]
        ru = rs[i]
        kl = np.sqrt(1j * omega * mu0 * conductivities[i-1])
        ku = np.sqrt(1j * omega * mu0 * conductivities[i+1])
        
        if np.abs(kj * rj) > n * cutoff_factor:
            juu, djuu, yuu, dyuu = get_spherical_harmonics(ku * ru, n, derivatives=True)
            del_l = -(juu + djuu / (1j * kj * rj)) / (yuu + dyuu / (1j * kj * rj))
            print('Layer {}: kr >> n'.format(i+1))
            if 1 / ku.imag > (ru - rl) / skin_depth_factor and 1 / ku.imag < (ru - rl) * skin_depth_factor:
                print('Approx. Invalid - 1/Im(k) = {}, layer thickness = {}'.format(1/ku.imag, ru-rl))
        
        elif np.abs(kj * rj) < n / cutoff_factor:
            jn1uu, yn1uu = get_spherical_harmonics(ku * ru, n + 1)
            jn_1uu, yn_1uu = get_spherical_harmonics(ku * ru, n - 1)
            jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
            jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
            AL = -(jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
            del_l = -(jn1uu - AL * jn_1uu * (rl / ru)**(2 * n + 1)) / (yn1uu - AL * yn_1uu * (rl / ru)**(2 * n + 1))
            print('Layer {}: kr << n'.format(i+1))
            if 1 / ku.imag > (ru - rl) / skin_depth_factor and 1 / ku.imag < (ru - rl) * skin_depth_factor:
                print('Approx. Invalid - 1/Im(k) = {}, layer thickness = {}'.format(1/ku.imag, ru-rl))

        else:
            jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
            jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

            beta_ul = jll * dyul - yul * djll
            gamma_ul = yll * dyul - yul * dyll
            delta_ul = jul * djll - jll * djul
            epsi_ul = jul * dyll - yll * djul
            
            del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)

    kJ = np.sqrt(1j * omega * mu0 * conductivities[-1])
    rJ = rs[-1]

    rl = rs[-2]

    if np.abs(kJ * rJ) > n * cutoff_factor:
        Ae = 1
        print('Layer {}: kr >> n'.format(i))
        if 1 / kJ.imag > (ru - rl) / skin_depth_factor and  1 / kJ.imag < (ru - rl) * skin_depth_factor:
            print('Approx. Invalid - 1/Im(k) = {}, layer thickness = {}'.format(1/ku.imag, ru-rl))

    elif np.abs(kJ * rJ) < n / cutoff_factor:
        kl = np.sqrt(1j * omega * mu0 * conductivities[-2])
        jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
        jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
        Ae = (rl / rJ)**(2 * n + 1) * (jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
        print('Layer {}: kr << n'.format(i))
        if 1 / kJ.imag > (ru - rl) / skin_depth_factor and  1 / kJ.imag < (ru - rl) * skin_depth_factor:
            print('Approx. Invalid - 1/Im(k) = {}, layer thickness = {}'.format(1/ku.imag, ru-rl))
    
    else:
        jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
        jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
        Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)

    return -Ae

def Aeiphi_Styczinski2(conductivities, rs, n, omega):
    def get_spherical_harmonics(kr, n, derivatives=False):
        jn = sps.spherical_jn(n, kr)
        yn = sps.spherical_yn(n, kr)
        if derivatives == True:
            d_jn = sps.spherical_jn(n, kr, derivative=True)
            d_yn = sps.spherical_yn(n, kr, derivative=True)
            return jn, d_jn, yn, d_yn
        return jn, yn

    # inner boundary
    k1 = np.sqrt(1j * omega * mu0 * conductivities[0])
    r1 = rs[0]
    k2 = np.sqrt(1j * omega * mu0 * conductivities[1])

    cutoff_factor = 10
    skin_depth_factor = 10

    if 1 / k1.imag > r1 / skin_depth_factor and 1 / k1.imag < r1 * skin_depth_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
        delta_ul = jul * djll - jll * djul
        beta_ul = jll * dyul - yul * djll
        del_l = delta_ul / beta_ul

    else:
        if np.abs(k1 * r1) > n * cutoff_factor:
            jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n, derivatives=True)
            del_l = -(jul + djul / (1j * k1 * r1)) / (yul + dyul / (1j * k1 * r1))
            print('Layer 1: kr >> n')    

        elif np.abs(k1 * r1) < n / cutoff_factor:
            jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
            del_l = -(jul / yul)
            print('Layer 1: kr << n')

        else:
            jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
            jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
            delta_ul = jul * djll - jll * djul
            beta_ul = jll * dyul - yul * djll
            del_l = delta_ul / beta_ul

    for i in range(1, len(rs) - 1):
        kj = np.sqrt(1j * omega * mu0 * conductivities[i])
        rj = rs[i]

        rl = rs[i-1]
        ru = rs[i]
        kl = np.sqrt(1j * omega * mu0 * conductivities[i-1])
        ku = np.sqrt(1j * omega * mu0 * conductivities[i+1])
        
        if 1 / ku.imag > (ru - rl) / skin_depth_factor and 1 / ku.imag < (ru - rl) * skin_depth_factor:
            jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
            jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

            beta_ul = jll * dyul - yul * djll
            gamma_ul = yll * dyul - yul * dyll
            delta_ul = jul * djll - jll * djul
            epsi_ul = jul * dyll - yll * djul
            
            del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)
        else:
            if np.abs(kj * rj) > n * cutoff_factor:
                juu, djuu, yuu, dyuu = get_spherical_harmonics(ku * ru, n, derivatives=True)
                del_l = -(juu + djuu / (1j * kj * rj)) / (yuu + dyuu / (1j * kj * rj))
                print('Layer {}: kr >> n'.format(i+1))

            elif np.abs(kj * rj) < n / cutoff_factor:
                jn1uu, yn1uu = get_spherical_harmonics(ku * ru, n + 1)
                jn_1uu, yn_1uu = get_spherical_harmonics(ku * ru, n - 1)
                jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
                jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
                AL = -(jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
                del_l = -(jn1uu - AL * jn_1uu * (rl / ru)**(2 * n + 1)) / (yn1uu - AL * yn_1uu * (rl / ru)**(2 * n + 1))
                print('Layer {}: kr << n'.format(i+1))

            else:
                jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
                jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

                beta_ul = jll * dyul - yul * djll
                gamma_ul = yll * dyul - yul * dyll
                delta_ul = jul * djll - jll * djul
                epsi_ul = jul * dyll - yll * djul
                
                del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)

    kJ = np.sqrt(1j * omega * mu0 * conductivities[-1])
    rJ = rs[-1]

    rl = rs[-2]

    if 1 / kJ.imag > (ru - rl) / skin_depth_factor and  1 / kJ.imag < (ru - rl) * skin_depth_factor:
        jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
        jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
        Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)
    
    else:
        if np.abs(kJ * rJ) > n * cutoff_factor:
            Ae = 1
            print('Layer {}: kr >> n'.format(i))

        elif np.abs(kJ * rJ) < n / cutoff_factor:
            kl = np.sqrt(1j * omega * mu0 * conductivities[-2])
            jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
            jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
            Ae = (rl / rJ)**(2 * n + 1) * (jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
            print('Layer {}: kr << n'.format(i))
        
        else:
            jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
            jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
            Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)

    return -Ae

def Aeiphi_Styczinski_many2(conductivities, rs, n, omega):
    def get_spherical_harmonics(kr, n, derivatives=False):
        jn = sps.spherical_jn(n, kr)
        yn = sps.spherical_yn(n, kr)
        if derivatives == True:
            d_jn = sps.spherical_jn(n, kr, derivative=True)
            d_yn = sps.spherical_yn(n, kr, derivative=True)
            return jn, d_jn, yn, d_yn
        return jn, yn

    # inner boundary
    k1 = np.sqrt(1j * omega * mu0 * conductivities[0])
    r1 = rs[0]
    k2 = np.sqrt(1j * omega * mu0 * conductivities[1])

    cutoff_factor = 10
    skin_depth_factor = 10

    if 1 / k1.imag > r1 / skin_depth_factor and 1 / k1.imag < r1 * skin_depth_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
        delta_ul = jul * djll - jll * djul
        beta_ul = jll * dyul - yul * djll
        del_l = delta_ul / beta_ul

    else:
        if np.abs(k1 * r1) > n * cutoff_factor:
            jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n, derivatives=True)
            del_l = -(jul + djul / (1j * k1 * r1)) / (yul + dyul / (1j * k1 * r1))  

        elif np.abs(k1 * r1) < n / cutoff_factor:
            jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
            del_l = -(jul / yul)
            
        else:
            jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
            jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
            delta_ul = jul * djll - jll * djul
            beta_ul = jll * dyul - yul * djll
            del_l = delta_ul / beta_ul

    for i in range(1, len(rs) - 1):
        kj = np.sqrt(1j * omega * mu0 * conductivities[i])
        rj = rs[i]

        rl = rs[i-1]
        ru = rs[i]
        kl = np.sqrt(1j * omega * mu0 * conductivities[i-1])
        ku = np.sqrt(1j * omega * mu0 * conductivities[i+1])
        
        if 1 / ku.imag > (ru - rl) / skin_depth_factor and 1 / ku.imag < (ru - rl) * skin_depth_factor:
            jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
            jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

            beta_ul = jll * dyul - yul * djll
            gamma_ul = yll * dyul - yul * dyll
            delta_ul = jul * djll - jll * djul
            epsi_ul = jul * dyll - yll * djul
            
            del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)
        else:
            if np.abs(kj * rj) > n * cutoff_factor:
                juu, djuu, yuu, dyuu = get_spherical_harmonics(ku * ru, n, derivatives=True)
                del_l = -(juu + djuu / (1j * kj * rj)) / (yuu + dyuu / (1j * kj * rj))

            elif np.abs(kj * rj) < n / cutoff_factor:
                jn1uu, yn1uu = get_spherical_harmonics(ku * ru, n + 1)
                jn_1uu, yn_1uu = get_spherical_harmonics(ku * ru, n - 1)
                jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
                jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
                AL = -(jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
                del_l = -(jn1uu - AL * jn_1uu * (rl / ru)**(2 * n + 1)) / (yn1uu - AL * yn_1uu * (rl / ru)**(2 * n + 1))

            else:
                jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
                jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

                beta_ul = jll * dyul - yul * djll
                gamma_ul = yll * dyul - yul * dyll
                delta_ul = jul * djll - jll * djul
                epsi_ul = jul * dyll - yll * djul
                
                del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)

    kJ = np.sqrt(1j * omega * mu0 * conductivities[-1])
    rJ = rs[-1]

    rl = rs[-2]

    if 1 / kJ.imag > (ru - rl) / skin_depth_factor and  1 / kJ.imag < (ru - rl) * skin_depth_factor:
        jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
        jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
        Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)
    
    else:
        if np.abs(kJ * rJ) > n * cutoff_factor:
            Ae = 1

        elif np.abs(kJ * rJ) < n / cutoff_factor:
            kl = np.sqrt(1j * omega * mu0 * conductivities[-2])
            jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
            jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
            Ae = (rl / rJ)**(2 * n + 1) * (jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
        
        else:
            jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
            jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
            Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)

    return -Ae

def Aeiphi_Styczinski_many(conductivities, rs, n, omega):
    def get_spherical_harmonics(kr, n, derivatives=False):
        jn = sps.spherical_jn(n, kr)
        yn = sps.spherical_yn(n, kr)
        if derivatives == True:
            d_jn = sps.spherical_jn(n, kr, derivative=True)
            d_yn = sps.spherical_yn(n, kr, derivative=True)
            return jn, d_jn, yn, d_yn
        return jn, yn

    # inner boundary
    k1 = np.sqrt(1j * omega * mu0 * conductivities[0])
    r1 = rs[0]
    k2 = np.sqrt(1j * omega * mu0 * conductivities[1])

    cutoff_factor = 100
    
    if np.abs(k1 * r1) > n * cutoff_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n, derivatives=True)
        del_l = -(jul + djul / (1j * k1 * r1)) / (yul + dyul / (1j * k1 * r1))

    elif np.abs(k1 * r1) < n / cutoff_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        del_l = -(jul / yul)

    else:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
        delta_ul = jul * djll - jll * djul
        beta_ul = jll * dyul - yul * djll
        del_l = delta_ul / beta_ul

    i = 2
    for i in range(1, len(rs) - 1):
        kj = np.sqrt(1j * omega * mu0 * conductivities[i])
        rj = rs[i]
        
        rl = rs[i-1]
        ru = rs[i]
        kl = np.sqrt(1j * omega * mu0 * conductivities[i-1])
        ku = np.sqrt(1j * omega * mu0 * conductivities[i+1])
        if np.abs(kj * rj) > n * cutoff_factor:
            juu, djuu, yuu, dyuu = get_spherical_harmonics(ku * ru, n)
            del_l = -(juu + djuu / (1j * kj * rj)) / (yuu + dyuu / (1j * kj * rj))
        
        elif np.abs(kj * rj) < n /cutoff_factor:
            jn1uu, yn1uu = get_spherical_harmonics(ku * ru, n + 1)
            jn_1uu, yn_1uu = get_spherical_harmonics(ku * ru, n - 1)
            jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
            jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
            AL = -(jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
            del_l = -(jn1uu - AL * jn_1uu * (rl / ru)**(2 * n + 1)) / (yn1uu - AL * yn_1uu * (rl / ru)**(2 * n + 1))

        else:
            jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
            jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

            beta_ul = jll * dyul - yul * djll
            gamma_ul = yll * dyul - yul * dyll
            delta_ul = jul * djll - jll * djul
            epsi_ul = jul * dyll - yll * djul
            
            del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)

    kJ = np.sqrt(1j * omega * mu0 * conductivities[-1])
    rJ = rs[-1]

    rl = rs[-2]

    if np.abs(kJ * rJ) > n * cutoff_factor:
        Ae = 1

    elif np.abs(kJ * rJ) < n / cutoff_factor:
        kl = np.sqrt(1j * omega * mu0 * conductivities[-2])
        jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
        jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
        Ae = (rl / rJ)**(2 * n + 1) * (jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
    
    else:
        jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
        jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
        Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)

    return -Ae
    
def Aeiphi_Styczinski_many3(conductivities, rs, n, omega):
    def get_spherical_harmonics(kr, n, derivatives=False):
        jn = sps.spherical_jn(n, kr)
        yn = sps.spherical_yn(n, kr)
        if derivatives == True:
            d_jn = sps.spherical_jn(n, kr, derivative=True)
            d_yn = sps.spherical_yn(n, kr, derivative=True)
            return jn, d_jn, yn, d_yn
        return jn, yn

    # inner boundary
    k1 = np.sqrt(1j * omega * mu0 * conductivities[0])
    r1 = rs[0]
    k2 = np.sqrt(1j * omega * mu0 * conductivities[1])

    cutoff_factor = 10
    skin_depth_factor = 10
    
    if np.abs(k1 * r1) > n * cutoff_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n, derivatives=True)
        del_l = -(jul + djul / (1j * k1 * r1)) / (yul + dyul / (1j * k1 * r1))

    elif np.abs(k1 * r1) < n / cutoff_factor:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        del_l = -(jul / yul)

    else:
        jul, djul, yul, dyul = get_spherical_harmonics(k2 * r1, n + 1, derivatives=True)
        jll, djll = get_spherical_harmonics(k1 * r1, n + 1)
        delta_ul = jul * djll - jll * djul
        beta_ul = jll * dyul - yul * djll
        del_l = delta_ul / beta_ul

    for i in range(1, len(rs) - 1):
        kj = np.sqrt(1j * omega * mu0 * conductivities[i])
        rj = rs[i]

        rl = rs[i-1]
        ru = rs[i]
        kl = np.sqrt(1j * omega * mu0 * conductivities[i-1])
        ku = np.sqrt(1j * omega * mu0 * conductivities[i+1])
        
        if np.abs(kj * rj) > n * cutoff_factor:
            juu, djuu, yuu, dyuu = get_spherical_harmonics(ku * ru, n, derivatives=True)
            del_l = -(juu + djuu / (1j * kj * rj)) / (yuu + dyuu / (1j * kj * rj))
        
        elif np.abs(kj * rj) < n / cutoff_factor:
            jn1uu, yn1uu = get_spherical_harmonics(ku * ru, n + 1)
            jn_1uu, yn_1uu = get_spherical_harmonics(ku * ru, n - 1)
            jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
            jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
            AL = -(jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
            del_l = -(jn1uu - AL * jn_1uu * (rl / ru)**(2 * n + 1)) / (yn1uu - AL * yn_1uu * (rl / ru)**(2 * n + 1))

        else:
            jul, djul, yul, dyul = get_spherical_harmonics(ku * rl, n + 1, derivatives=True)
            jll, djll, yll, dyll = get_spherical_harmonics(kl * rl, n + 1, derivatives=True)

            beta_ul = jll * dyul - yul * djll
            gamma_ul = yll * dyul - yul * dyll
            delta_ul = jul * djll - jll * djul
            epsi_ul = jul * dyll - yll * djul
            
            del_l = (delta_ul + del_l * epsi_ul) / (beta_ul + del_l * gamma_ul)

    kJ = np.sqrt(1j * omega * mu0 * conductivities[-1])
    rJ = rs[-1]

    rl = rs[-2]

    if np.abs(kJ * rJ) > n * cutoff_factor:
        Ae = 1

    elif np.abs(kJ * rJ) < n / cutoff_factor:
        kl = np.sqrt(1j * omega * mu0 * conductivities[-2])
        jn1ll, yn1ll = get_spherical_harmonics(kl * rl, n + 1)
        jn_1ll, yn_1ll = get_spherical_harmonics(kl * rl, n - 1)
        Ae = (rl / rJ)**(2 * n + 1) * (jn1ll + del_l * yn1ll) / (jn_1ll + del_l * yn_1ll)
    
    else:
        jn1R, yn1R = get_spherical_harmonics(kJ * rJ, n + 1)
        jn_1R, yn_1R = get_spherical_harmonics(kJ * rJ, n - 1)
        Ae = -(jn1R + del_l * yn1R) / (jn_1R + del_l * yn_1R)

    return -Ae

def B_induced_aeiphi_minimiser(O, B_external, t_longperiod, omega, A):
    """
    Calculate the induced magnetic field with finite conductivity
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param omega: angular frequency of inducing field in
    :param conductivities: array of conductivities of the layers in S
    :param radii: array of radii of the layers in m
    :return: time evolution array of Bx, By, Bz in nT
    """
    orbit = O.copy()
    t_red = orbit[0]
    orbit = orbit.transpose()
    B_ext = B_external

    
    print('Aeiphi = {}'.format(A))

    t = t_longperiod
    
    # print(A)
    phi = -np.angle(A)
    # print(phi)
    # print(omega)
    t_phi = phi / omega
    # print(t_phi)
    t_interval = t[1] - t[0]
    # print(t_interval)
    n_intervals = int(np.round(t_phi / t_interval))
    B_ext = np.roll(B_ext, n_intervals, axis=0)

    tB_ext = np.c_[t_longperiod, B_ext]
    tB_ext = tB_ext.transpose()

    tB_reduced = []
    for i in range(len(t_red)):
        index = find_nearest_index(tB_ext[0], t_red[i])
        tB_i = tB_ext[:, int(index)]
        tB_reduced.append(tB_i)
    tB_reduced = np.transpose(tB_reduced)
    B_reduced = tB_reduced[1:]
    B_ext = B_reduced.transpose()

    # _M = M = -(2 * pi / mu0) * A * (radii[-1]**3)
    _M = M = -(2 * pi / mu0) * A * (R_C**3)
    # print('M = {}'.format(_M))

    Bind_evolution = []
    for B_ext, vector in zip(B_ext, orbit):

        position = vector[1:4]
        # M = -(2 * pi / mu0) * A * B_ext * (radii[-1]**3)
        M = B_ext * _M

        rmag = np.linalg.norm(position)
        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind = Bind.real
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)

        

        
    
