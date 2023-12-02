import numpy as np
n = 1


def B_sph_components(r,
                     theta,
                     phi,
                     coeffs,
                     m_ints,
                     leg_poly: dict,
                     dleg_poly: dict):
    Br = 0
    Btheta = 0
    Bphi = 0
    for m in m_ints:
        if m <= n:
            num_r = (n + 1) * ((1 / r) ** (n + 2))
            num_theta = (1 / r) ** (n + 2)
            num_phi = m * (1 / r) ** (n + 2)

            lp = leg_poly[f"P{n}{m}"]
            dlp = dleg_poly[f"dP{n}{m}"]

            gc_r_theta = coeffs[f"g{n}{m}"] * np.cos(m * phi) + coeffs[f"h{n}{m}"] * np.sin(m * phi)
            gc_phi = coeffs[f"g{n}{m}"] * np.sin(m * phi) - coeffs[f"h{n}{m}"] * np.cos(m * phi)

            Br += (num_r * lp * gc_r_theta)
            Btheta += (num_theta * dlp * gc_r_theta)
            Bphi += (num_phi * lp * gc_phi)

    Btheta = - Btheta
    Bphi = 1 / (np.sin(theta)) * Bphi

    return Br, Btheta, Bphi