from math import *
import numpy as np


class Galcoords:

    # Initializer / Instance Attributes
    def __init__(self, usun=14., vrot_sun=253., wsun=7., rsun=8.122, **kwargs):
        self.usun = usun
        self.vrot_sun = vrot_sun
        self.wsun = wsun
        self.rsun = rsun

        self.lb = False
        self.uvw = False
        self.gal_pos = False
        self.gal_vel = False
        self.sgr_LB = False
        self.sgr_pos = False
        self.sgr_vel = False

        self.ra = None
        self.dec = None
        self.l = None
        self.b = None
        self.pm_ra = None
        self.pm_ra_e = None
        self.pm_dec = None
        self.pm_dec_e = None
        self.v_rad = None
        self.v_rad_e = None
        self.dist = None
        self.sigma_dist = None
        self.parallax = None
        self.parallax_e = None

        self.update_values(**kwargs)

        if self.l is not None and self.b is not None:
            self.lb = True

    def update_values(self, ra=None, dec=None, l=None, b=None, pm_ra=None,
                      pm_ra_e=None, pm_dec=None, pm_dec_e=None, v_rad=None,
                      v_rad_e=None, dist=None, sigma_dist=None, parallax=None,
                      parallax_e=None, use_dist=False):
        '''
        parallaxes and proper motions should be entered in arcsec and arcsec/yr
        distances should be entered in kpc.  If distances and parallaxes are
        entered, parallaxes will be preferred and if lower and upper
        uncertainties are given on distances in addition to a mean uncertainty
        the distance uncertainty will be calculated using the average of the
        upper and lower.
        '''

        keys = ['ra', 'dec', 'l', 'b', 'pm_ra', 'pm_ra_e', 'pm_dec', 'pm_dec_e',
                'v_rad', 'v_rad_e', 'dist', 'sigma_dist', 'parallax',
                'parallax_e']

        input_values = [ra, dec, l, b, pm_ra, pm_ra_e, pm_dec, pm_dec_e, v_rad,
                        v_rad_e, dist, sigma_dist, parallax, parallax_e]

        for key, value in zip(keys, input_values):
            if value is not None:
                setattr(self, key, value)

        if use_dist:
            self.parallax = 1. / self.dist
            self.parallax_e = self.sigma_dist / (self.dist**2.)

    def check_inputs(self, input_param, input_values, uncertainty_keys,
                     input_uncertainties):

        for key, value in input_param.items():
            if value is None:
                raise TypeError(f'Please input a value for {key}.')

        for input_uncertainty, unc_key, input_value in zip(input_uncertainties,
                                                           uncertainty_keys,
                                                           input_values):
            if input_uncertainty is None:
                setattr(self, unc_key, np.zeros(len(input_value)))

    def calculate_helio_vel(self, **kwargs):

        self.update_values(**kwargs)

        input_param = {'ra': self.ra, 'dec': self.dec,
                       'parallax or dist': self.parallax,
                       'v_rad': self.v_rad, 'pm_ra': self.pm_ra,
                       'pm_dec': self.pm_dec}
        input_values = [self.parallax, self.v_rad, self.pm_ra, self.pm_dec]
        uncertainty_keys = ['parallax_e', 'v_rad_e', 'pm_ra_e', 'pm_dec_e']
        input_uncertainties = [self.parallax_e, self.v_rad_e, self.pm_ra_e,
                               self.pm_dec_e]

        self.check_inputs(input_param, input_values, uncertainty_keys,
                          input_uncertainties)

        self._calc_galuvw()

    def calculate_gal_lb(self, **kwargs):

        self.update_values(**kwargs)
        input_param = {'ra': self.ra, 'dec': self.dec}
        input_values = []
        uncertainty_keys = []
        input_uncertainties = []
        self.check_inputs(input_param, input_values, uncertainty_keys,
                          input_uncertainties)
        self._calc_lb()
        self.lb = True

    def calculate_gal_pos(self, **kwargs):

        if self.lb:
            self.update_values(**kwargs)
        else:
            self.calculate_gal_lb(**kwargs)

        input_param = {'dist': self.dist}
        input_values = [self.dist]
        uncertainty_keys = ['sigma_dist']
        input_uncertainties = [self.sigma_dist]
        self.check_inputs(input_param, input_values, uncertainty_keys,
                          input_uncertainties)
        self._calc_galpos()
        self.gal_pos = True

    def calculate_gal_vel(self, update_inputs=False, **kwargs):

        if not self.uvw or update_inputs:
            self.calculate_helio_vel(**kwargs)
        if not self.gal_pos or update_inputs:
            self.calculate_gal_pos(**kwargs)

        self._calc_galvel()
        self.gal_vel = True

    def calculate_sgr_lb(self, **kwargs):

        if self.gal_pos:
            self._calc_sgrlb_xyz()
        else:
            if not self.lb:
                self.calculate_gal_lb(**kwargs)
            self._calc_sgrlb_lb()
        self.sgr_LB = True

    def calculate_sgr_pos(self, **kwargs):

        if not self.sgr_LB:
            self.calculate_sgr_lb(**kwargs)
        if not self.gal_pos:
            self.calculate_gal_pos(**kwargs)

        self._calc_sgrpos()
        self.sgr_pos = True

    def calculate_sgr_system(self, **kwargs):

        if not self.sgr_LB:
            self.calculate_sgr_lb(**kwargs)
        if not self.gal_vel:
            self.calculate_gal_vel(**kwargs)

        self._calc_sgrvel()
        self.sgr_vel = True

    def _calc_galuvw(self):

        self.dist, self.sigma_dist, self.U, self.sigma_U, self.V, \
            self.sigma_V, self.W, self.sigma_W = gal_uvw(self.ra, self.dec,
                                                         self.parallax,
                                                         self.parallax_e,
                                                         self.v_rad,
                                                         self.v_rad_e,
                                                         self.pm_ra,
                                                         self.pm_ra_e,
                                                         self.pm_dec,
                                                         self.pm_dec_e)

    def _calc_lb(self):

        self.l, self.b = gal_lb(self.ra, self.dec)

    def _calc_galpos(self):

        self.x, self.y, self.z, self.R_cyn, self.phi, self.r_sph, \
            self.sigma_x, self.sigma_y, self.sigma_z, self.sigma_R_cyn, \
            self.sigma_r_sph, self.sigma_phi = gal_coords_err(self.l, self.b,
                                                              self.dist,
                                                              self.sigma_dist,
                                                              self.rsun)

    def _calc_galvel(self):

        self.vx, self.sigma_vx, self.vy, self.sigma_vy, self.vz, \
            self.sigma_vz, self.vR, self.sigma_vR, self.vphi, \
            self.sigma_vphi = gal_vel(self.R_cyn, self.phi, self.sigma_phi,
                                      self.usun, self.vrot_sun, self.wsun,
                                      self.U, self.sigma_U, self.V,
                                      self.sigma_V, self.W, self.sigma_W)

    def _calc_sgrlb_xyz(self):

        self.lambda_sun, self.beta_sun = sgr_coords(self.x, self.y, self.z,
                                                    self.rsun)

    def _calc_sgrlb_lb(self):

        self.lambda_sun, self.beta_sun = sgr_coords_lb(self.l, self.b)

    def _calc_sgrpos(self):
        self.lambda_gc, self.beta_gc, self.xs, self.ys, self.zs, self.r_cys, \
            self.sigma_lambda_gc, self.sigma_xs, self.sigma_ys, self.sigma_zs, \
            self.sigma_r_cys = sgr_pos(self.x, self.y, self.z, self.rsun,
                                       self.sigma_x, self.sigma_y, self.sigma_z)

    def _calc_sgrvel(self):

        self.lambda_gc, self.beta_gc, self.xs, self.ys, self.zs, self.r_cys, \
            self.vxs, self.vys, self.vzs, self.vrs, self.vphis, \
            self.sigma_lambda_gc, self.sigma_xs, self.sigma_ys, self.sigma_zs, \
            self.sigma_r_cys, self.sigma_vxs, self.sigma_vys, self.sigma_vzs, \
            self.sigma_vrs, self.sigma_vphis = sgr_system(self.x, self.y,
                                                          self.z, self.vx,
                                                          self.vy, self.vz,
                                                          self.rsun,
                                                          self.sigma_x,
                                                          self.sigma_y,
                                                          self.sigma_z,
                                                          self.sigma_vx,
                                                          self.sigma_vy,
                                                          self.sigma_vz)


def gal_uvw(ra, dec, parallax, parallax_e, v_r, v_r_e, pm_ra, pm_ra_e,
            pm_dec, pm_dec_e):
    '''
    gal_uvw calculates distances, and U, V, and W velocities and their
    uncertainties from sky positions (RA and Dec), parallaxes, radial
    velocities, and proper motions following calculations from
    Johnson and Soderblom (1987, AJ, 93, 864)

    Requires that ra, dec, parallax, radial velocity (helio/barycentric; v_r)
    proper motion in ra and dec (pm_ra and pm_dec respectively) and their
    uncertainties (i.e., variables ending in _e) are given.  If uncertainties
    are not known, enter zeros.

    Entries can be individaul values or numpy-like arrays of matching lengths.
    '''

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


def sgr_pos(x, y, z, xsun, sigma_x, sigma_y, sigma_z):

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

    return lambda_s, beta, xs, ys, zs, r_cys, lambda_s_err, xs_err, ys_err, \
        zs_err, r_cys_err


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
