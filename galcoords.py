from math import *
import numpy as np


def gal_uvw(ra, dec, parallax, parallax_e, v_r, v_r_e, pm_ra, pm_ra_e,
            pm_dec, pm_dec_e):
    # theta_0 = np.deg2rad(123.) #B1950
    # ra_ngp = np.deg2rad(192.25) #B1950
    # dec_ngp = np.deg2rad(27.4) #B1950
    theta_0 = np.deg2rad(122.93)  # J2000
    ra_ngp = np.deg2rad(192.85)  # J2000
    dec_ngp = np.deg2rad(27.133333333)  # J2000

    k = 4.74057

    t1 = np.zeros((3, 3))
    t2 = np.zeros((3, 3))
    t3 = np.zeros((3, 3))
    a1 = np.zeros((3, 3))
    a2 = np.zeros((3, 3))

    t1[0, 0] = np.cos(theta_0)
    t1[0, 1] = np.sin(theta_0)
    t1[1, 0] = np.sin(theta_0)
    t1[1, 1] = -1. * np.cos(theta_0)
    t1[2, 2] = 1.

    t2[0, 0] = -1. * np.sin(dec_ngp)
    t2[0, 2] = np.cos(dec_ngp)
    t2[2, 0] = np.cos(dec_ngp)
    t2[1, 1] = -1.
    t2[2, 2] = np.sin(dec_ngp)

    t3[0, 0] = np.cos(ra_ngp)
    t3[0, 1] = np.sin(ra_ngp)
    t3[1, 0] = np.sin(ra_ngp)
    t3[1, 1] = -1. * np.cos(ra_ngp)
    t3[2, 2] = 1.

    T = np.dot(np.dot(t1, t2), t3)

    t00 = T[0, 0]
    t01 = T[0, 1]
    t02 = T[0, 2]
    t10 = T[1, 0]
    t11 = T[1, 1]
    t12 = T[1, 2]
    t20 = T[2, 0]
    t21 = T[2, 1]
    t22 = T[2, 2]

    cosa = np.cos(np.radians(ra))
    sina = np.sin(np.radians(ra))
    cosd = np.cos(np.radians(dec))
    sind = np.sin(np.radians(dec))

    a00 = cosa * cosd
    a01 = -1. * sina
    a02 = -1. * cosa * sind
    a10 = sina * cosd
    a11 = 1 * cosa
    a12 = -1. * sina * sind
    a20 = 1. * sind
    a21 = 0.
    a22 = 1. * cosd

    b00 = t00 * a00 + t01 * a10 + t02 * a20
    b01 = t00 * a01 + t01 * a11 + t02 * a21
    b02 = t00 * a02 + t01 * a12 + t02 * a22
    b10 = t10 * a00 + t11 * a10 + t12 * a20
    b11 = t10 * a01 + t11 * a11 + t12 * a21
    b12 = t10 * a02 + t11 * a12 + t12 * a22
    b20 = t20 * a00 + t21 * a10 + t22 * a20
    b21 = t20 * a01 + t21 * a11 + t22 * a21
    b22 = t20 * a02 + t21 * a12 + t22 * a22

    # Calculate UVW

    v_a = k * pm_ra / parallax
    v_d = k * pm_dec / parallax

    U = b00 * v_r + b01 * v_a + b02 * v_d
    V = b10 * v_r + b11 * v_a + b12 * v_d
    W = b20 * v_r + b21 * v_a + b22 * v_d

    # Calculate sigma squared uvw

    c00 = b00**2.
    c01 = b01**2.
    c02 = b02**2.
    c10 = b10**2.
    c11 = b11**2.
    c12 = b12**2.
    c20 = b20**2.
    c21 = b21**2.
    c22 = b22**2.

    var_vr = v_r_e**2.
    var_va = (k / parallax)**2. * (pm_ra_e**2. +
                                   (pm_ra * parallax_e / parallax)**2.)
    var_vd = (k / parallax)**2. * (pm_dec_e**2. +
                                   (pm_dec * parallax_e / parallax)**2.)

    sigma_U2 = (c00 * var_vr + c01 * var_va + c02 * var_vd) + \
               (2 * pm_ra * pm_dec * k**2. * parallax_e**2. / parallax**4.) * \
               (b01 * b02)
    sigma_V2 = (c10 * var_vr + c11 * var_va + c12 * var_vd) + \
               (2 * pm_ra * pm_dec * k**2. * parallax_e**2. / parallax**4.) * \
               (b11 * b12)
    sigma_W2 = (c20 * var_vr + c21 * var_va + c22 * var_vd) + \
               (2 * pm_ra * pm_dec * k**2. * parallax_e**2. / parallax**4.) * \
               (b21 * b22)

    sigma_U = np.sqrt(sigma_U2)
    sigma_V = np.sqrt(sigma_V2)
    sigma_W = np.sqrt(sigma_W2)

    dist = 1. / parallax
    sigma_dist = parallax_e / parallax**2.
    dist_kpc = dist / 1000.
    sigma_dist_kpc = sigma_dist / 1000.

    return dist_kpc, sigma_dist_kpc, U, sigma_U, V, sigma_V, W, sigma_W


def gal_coords(l, b, dist, rsun):
    l_rad = np.radians(l)
    b_rad = np.radians(b)

    z = dist * np.sin(b_rad)
    y = dist * np.cos(b_rad) * np.sin(l_rad)
    x = rsun - dist * np.cos(b_rad) * np.cos(l_rad)

    R_cyn = np.sqrt(x**2. + y**2.)
    r_sph = np.sqrt(x**2. + y**2. + z**2.)
    phi = np.arctan2(y, x)

    return x, y, z, R_cyn, phi, r_sph


def gal_coords_err(l, b, dist, dist_err, rsun):
    l_rad = np.radians(l)
    b_rad = np.radians(b)

    z = dist * np.sin(b_rad)
    y = dist * np.cos(b_rad) * np.sin(l_rad)
    x = rsun - dist * np.cos(b_rad) * np.cos(l_rad)

    z_err = dist_err * np.sin(b_rad)
    y_err = dist_err * np.cos(b_rad) * np.sin(l_rad)
    x_err = dist_err * np.cos(b_rad) * np.cos(l_rad)

    R_cyn = np.sqrt(x**2. + y**2.)
    r_sph = np.sqrt(x**2. + y**2. + z**2.)
    phi = np.arctan2(y, x)

    R_cyn_err = np.sqrt((x_err**2. * x**2. + y_err**2. * y**2.) / R_cyn**2.)
    r_sph_err = np.sqrt((x_err**2. * x**2. + y_err**2. * y **
                         2. + z_err**2. * z**2.) / r_sph**2.)
    phi_err = np.sqrt((x_err**2. * y**2. + y_err**2. * x**2.) / R_cyn**4.)

    return x, y, z, R_cyn, phi, r_sph, x_err, y_err, z_err, R_cyn_err, \
        r_sph_err, phi_err


def gal_vel(R_cyn, phi, sigma_phi, Usun, Vrot_sun, Wsun, U, sigma_U, V,
            sigma_V, W, sigma_W):
    # Vx, U and vR are measured as increasing toward the galactic center.
    # the anticenter is in the direction of phi = 0 or increasing x

    vx = U + Usun
    vy = V + Vrot_sun
    vz = W + Wsun

    sigma_vx = sigma_U
    sigma_vy = sigma_V
    sigma_vz = sigma_W

    vR = vx * np.cos(phi) - vy * np.sin(phi)
    vphi = vx * np.sin(phi) + vy * np.cos(phi)

    sigma_vR = np.sqrt((sigma_vx * np.cos(phi))**2. + (sigma_vy * np.sin(phi))
                       ** 2. + sigma_phi**2. * (-vx * np.sin(phi) - vy *
                                                np.cos(phi))**2.)
    sigma_vphi = np.sqrt((sigma_vx * np.sin(phi))**2. + (sigma_vy * np.cos(phi))
                         ** 2. + sigma_phi**2. * (vx * np.cos(phi) - vy *
                                                  np.sin(phi))**2.)

    return vx, sigma_vx, vy, sigma_vy, vz, sigma_vz, vR, sigma_vR, vphi, \
        sigma_vphi


def sgr_system(x, y, z, vx, vy, vz, xsun, sigma_x, sigma_y, sigma_z, sigma_vx,
               sigma_vy, sigma_vz):

    phi = np.deg2rad(180. + 3.75)
    theta = np.deg2rad(90. - 13.46)
    psiGC = np.deg2rad(180. + 21.604339)
    psi_sun = np.deg2rad(180. + 21.604339)

    # ang is the rotation of phiGC past 180
    ang = np.deg2rad(21.604399)

    xcenter = -8.5227
    ycenter = -0.3460
    zcenter = -0.828

    GCrot11 = np.cos(psiGC) * np.cos(phi) - np.cos(theta) * \
        np.sin(phi) * np.sin(psiGC)
    GCrot12 = np.cos(psiGC) * np.sin(phi) + np.cos(theta) * \
        np.cos(phi) * np.sin(psiGC)
    GCrot13 = np.sin(psiGC) * np.sin(theta)
    GCrot21 = -np.sin(psiGC) * np.cos(phi) - np.cos(theta) * \
        np.sin(phi) * np.cos(psiGC)
    GCrot22 = -np.sin(psiGC) * np.sin(phi) + np.cos(theta) * \
        np.cos(phi) * np.cos(psiGC)
    GCrot23 = np.cos(psiGC) * np.sin(theta)
    GCrot31 = np.sin(theta) * np.sin(phi)
    GCrot32 = -np.sin(theta) * np.cos(phi)
    GCrot33 = np.cos(theta)

    x = -x
    x = x + xsun

    temp = GCrot11 * (x + xcenter) + GCrot12 * \
        (y - ycenter) + GCrot13 * (z - zcenter)
    temp2 = GCrot21 * (x + xcenter) + GCrot22 * \
        (y - ycenter) + GCrot23 * (z - zcenter)
    zs = GCrot31 * (x + xcenter) + GCrot32 * \
        (y - ycenter) + GCrot33 * (z - zcenter)
    d = np.sqrt(temp * temp + temp2 * temp2 + zs * zs)

    temp_err = np.sqrt((GCrot11 * sigma_x)**2. +
                       (GCrot12 * sigma_y)**2. + (GCrot13 * sigma_z)**2.)
    temp2_err = np.sqrt((GCrot21 * sigma_x)**2. +
                        (GCrot22 * sigma_y)**2. + (GCrot23 * sigma_z)**2.)
    zs_err = np.sqrt((GCrot31 * sigma_x)**2. +
                     (GCrot32 * sigma_y)**2. + (GCrot33 * sigma_z)**2.)

    zs = -zs

    temp3 = np.rad2deg(np.arctan2(temp2, temp))
    temp3[temp3 < 0.] = temp3[temp3 < 0.] + 360.
    temp3 = temp3 + np.rad2deg(ang)
    temp3[temp3 > 360.] = temp3[temp3 > 360.] - 360.
    lambda_s = temp3
    beta = np.rad2deg(
        np.arcsin(zs / np.sqrt(temp * temp + temp2 * temp2 + zs * zs)))

    lambda_s_err = np.sqrt((temp_err**2. * temp2**2. +
                            temp2_err**2. * temp**2.) /
                           (temp**2. + temp2**2.)**2.)

    xs = temp * np.cos(ang) - temp2 * np.sin(ang)
    ys = temp * np.sin(ang) + temp2 * np.cos(ang)
    r_cys = np.sqrt(xs**2. + ys**2.)

    xs_err = np.sqrt((temp_err * np.cos(ang))**2. +
                     (temp2_err * np.sin(ang))**2.)
    ys_err = np.sqrt((temp_err * np.sin(ang))**2. +
                     (temp2_err * np.cos(ang))**2.)
    r_cys_err = np.sqrt(((xs * xs_err)**2. + (ys * ys_err)**2.) / r_cys**2.)

    vtemp = GCrot11 * vx + GCrot12 * vy + GCrot13 * vz
    vtemp2 = GCrot21 * vx + GCrot22 * vy + GCrot23 * vz
    vzs = GCrot31 * vx + GCrot32 * vy + GCrot33 * vz

    vtemp_err = np.sqrt((GCrot11 * sigma_vx)**2. +
                        (GCrot12 * sigma_vy)**2. + (GCrot13 * sigma_vz)**2.)
    vtemp2_err = np.sqrt((GCrot21 * sigma_vx)**2. +
                         (GCrot22 * sigma_vy)**2. + (GCrot23 * sigma_vz)**2.)
    vzs_err = np.sqrt((GCrot31 * sigma_vx)**2. +
                      (GCrot32 * sigma_vy)**2. + (GCrot33 * sigma_vz)**2.)

    vzs = -vzs
    vxs = vtemp * np.cos(ang) - vtemp2 * np.sin(ang)
    vys = vtemp * np.sin(ang) + vtemp2 * np.cos(ang)

    vxs_err = np.sqrt((vtemp_err * np.cos(ang))**2. +
                      (vtemp2_err * np.sin(ang))**2.)
    vys_err = np.sqrt((vtemp_err * np.sin(ang))**2. +
                      (vtemp2_err * np.cos(ang))**2.)

    vrs = vxs * np.cos(np.deg2rad(lambda_s)) + vys * \
        np.sin(np.deg2rad(lambda_s))
    vphis = vxs * np.sin(np.deg2rad(lambda_s)) - vys * \
        np.cos(np.deg2rad(lambda_s))

    sigma_vrs = np.sqrt((vxs_err * np.cos(np.radians(lambda_s)))**2. +
                        (vys_err * np.sin(np.radians(lambda_s))) ** 2. +
                        lambda_s_err**2. *
                        (-vxs * np.sin(np.radians(lambda_s)) -
                         vys * np.cos(np.radians(lambda_s)))**2.)
    sigma_vphis = np.sqrt((vxs_err * np.sin(np.radians(lambda_s)))**2. +
                          (vys_err * np.cos(np.radians(lambda_s))) ** 2. +
                          lambda_s_err**2. *
                          (vxs * np.cos(np.radians(lambda_s)) -
                           vys * np.sin(np.radians(lambda_s)))**2.)

    lambda_s_err = np.degrees(lambda_s_err)

    return lambda_s, beta, xs, ys, zs, r_cys, vxs, vys, vzs, vrs, vphis, \
        lambda_s_err, xs_err, ys_err, zs_err, r_cys_err, vxs_err, vys_err, \
        vzs_err, sigma_vrs, sigma_vphis


def sgr_coords_lb(l, b):

    x, y, z, r1, r2, phi_g = gal_coords(l, b, np.ones(len(l)), 0.)

    phi = np.deg2rad(180. + 3.75)
    theta = np.deg2rad(90. - 13.46)
    psi_s = np.deg2rad(180. + 14.111534)

    rot11 = np.cos(psi_s) * np.cos(phi) - np.cos(theta) * \
        np.sin(phi) * np.sin(psi_s)
    rot12 = np.cos(psi_s) * np.sin(phi) + np.cos(theta) * \
        np.cos(phi) * np.sin(psi_s)
    rot13 = np.sin(psi_s) * np.sin(theta)
    rot21 = -np.sin(psi_s) * np.cos(phi) - np.cos(theta) * \
        np.sin(phi) * np.cos(psi_s)
    rot22 = -np.sin(psi_s) * np.sin(phi) + np.cos(theta) * \
        np.cos(phi) * np.cos(psi_s)
    rot23 = np.cos(psi_s) * np.sin(theta)
    rot31 = np.sin(theta) * np.sin(phi)
    rot32 = -np.sin(theta) * np.cos(phi)
    rot33 = np.cos(theta)

    x = -x

    Xs = rot11 * x + rot12 * y + rot13 * z
    Ys = rot21 * x + rot22 * y + rot23 * z
    Zs = rot31 * x + rot32 * y + rot33 * z
    r = np.sqrt(Xs * Xs + Ys * Ys + Zs * Zs)

    Zs = -Zs

    lambda_sun = np.rad2deg(np.arctan2(Ys, Xs))
    lambda_sun[lambda_sun < 0.] = lambda_sun[lambda_sun < 0.] + 360.
    beta = np.rad2deg(np.arcsin(Zs / np.sqrt(Xs * Xs + Ys * Ys + Zs * Zs)))
    return lambda_sun, beta


def sgr_coords(x, y, z, xsun):

    phi = np.deg2rad(180. + 3.75)
    theta = np.deg2rad(90. - 13.46)
    psi_s = np.deg2rad(180. + 14.111534)

    rot11 = np.cos(psi_s) * np.cos(phi) - np.cos(theta) * \
        np.sin(phi) * np.sin(psi_s)
    rot12 = np.cos(psi_s) * np.sin(phi) + np.cos(theta) * \
        np.cos(phi) * np.sin(psi_s)
    rot13 = np.sin(psi_s) * np.sin(theta)
    rot21 = -np.sin(psi_s) * np.cos(phi) - np.cos(theta) * \
        np.sin(phi) * np.cos(psi_s)
    rot22 = -np.sin(psi_s) * np.sin(phi) + np.cos(theta) * \
        np.cos(phi) * np.cos(psi_s)
    rot23 = np.cos(psi_s) * np.sin(theta)
    rot31 = np.sin(theta) * np.sin(phi)
    rot32 = -np.sin(theta) * np.cos(phi)
    rot33 = np.cos(theta)

    x = -x
    x = x + xsun

    Xs = rot11 * x + rot12 * y + rot13 * z
    Ys = rot21 * x + rot22 * y + rot23 * z
    Zs = rot31 * x + rot32 * y + rot33 * z
    r = np.sqrt(Xs * Xs + Ys * Ys + Zs * Zs)

    Zs = -Zs

    lambda_sun = np.rad2deg(np.arctan2(Ys, Xs))
    lambda_sun[lambda_sun < 0.] = lambda_sun[lambda_sun < 0.] + 360.
    beta = np.rad2deg(np.arcsin(Zs / np.sqrt(Xs * Xs + Ys * Ys + Zs * Zs)))
    return lambda_sun, beta


def gal_lb(ra, dec):

    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    a0 = np.deg2rad(282.85)
    angp = np.deg2rad(192.8595)
    dngp = np.deg2rad(27.1283)
    agc = np.deg2rad(266.4050)
    dgc = np.deg2rad(-28.9362)

    #l0 = np.deg2rad(90+32.93)
    l0 = np.arccos(np.sin(dgc) * np.cos(dngp) - np.cos(dgc)
                   * np.sin(dngp) * np.cos(angp - agc))

    b = np.rad2deg(np.arcsin(np.sin(dec_rad) * np.sin(dngp) +
                             np.cos(dec_rad) * np.cos(dngp) *
                             np.cos(ra_rad - angp)))
    l = np.degrees(l0 - np.arctan2((np.sin(ra_rad - angp) * np.cos(dec_rad)),
                                   (np.sin(dec_rad) * np.cos(dngp) -
                                    np.cos(dec_rad) * np.sin(dngp) *
                                    np.cos(ra_rad - angp))))

    l[l < 0.] = l[l < 0.] + 360.

    return l, b
