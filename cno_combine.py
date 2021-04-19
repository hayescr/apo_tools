import numpy as np


def cno_combine(fe_h, c_fe, n_fe, o_fe, fe_h_err, c_fe_err, n_fe_err, o_fe_err):

    c_h = c_fe + fe_h
    n_h = n_fe + fe_h
    o_h = o_fe + fe_h

    # CALCULATE C+N and C+N+O VALUES
    # Using GREVESSE 2007/Asplund 2005 values for C, N, O, and Fe
    ch_solar = 8.39
    nh_solar = 7.78
    oh_solar = 8.66
    cnh_solar = np.log10(10.**ch_solar + 10.**nh_solar)
    cnoh_solar = np.log10(10.**ch_solar + 10.**nh_solar + 10.**oh_solar)

    ch = c_h + ch_solar
    nh = n_h + nh_solar
    oh = o_h + oh_solar

    # [C+N/Fe]
    cnh = np.log10(10.**ch + 10.**nh) - cnh_solar
    cn_fe = cnh - fe_h

    # [C+N+O/Fe]
    cnoh = np.log10(10.**ch + 10.**nh + 10.**oh) - cnoh_solar
    cno_fe = cnoh - fe_h

    # Errors

    c_h_err = np.sqrt(c_fe_err**2. + fe_h_err**2.)
    n_h_err = np.sqrt(n_fe_err**2. + fe_h_err**2.)
    o_h_err = np.sqrt(o_fe_err**2. + fe_h_err**2.)

    cn_ddch = 10.**ch / (10.**ch + 10.**nh)
    cn_ddnh = 10.**nh / (10.**ch + 10.**nh)

    cno_ddch = 10.**ch / (10.**ch + 10.**nh + 10.**oh)
    cno_ddnh = 10.**nh / (10.**ch + 10.**nh + 10.**oh)
    cno_ddoh = 10.**oh / (10.**ch + 10.**nh + 10.**oh)

    cn_fe_err = np.sqrt((cn_ddch * c_h_err)**2 +
                        (cn_ddnh * n_h_err)**2 + c_fe_err**2.)

    cno_fe_err = np.sqrt((cno_ddch * c_h_err)**2 + (cno_ddnh * n_h_err)
                         ** 2 + (cno_ddoh * o_h_err)**2 + c_fe_err**2.)

    return cn_fe, cno_fe, cn_fe_err, cno_fe_err
