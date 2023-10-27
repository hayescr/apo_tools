import numpy as np
from pkg_resources import resource_filename
filepath = resource_filename('apo_tools', 'prantzos2020_isotope_data.txt')

def calculate_isotope_fractions(z, isotopes, isotope_abund):
    elements = np.unique(z)
    isotope_dict = {}
    for element in elements:
        selection = z==element
        total_mass = isotope_abund[selection].sum()
        isotope_dict.update({element+iso/1000:np.round(abund/total_mass,5) for iso, abund in zip(isotopes[selection], isotope_abund[selection])})

    return isotope_dict

def create_iso_dict(rprocess=1., sprocess=1.):
    '''
    Calculates the isotopic ratios of custom r-/s-process mix.  Input values are set to
        scale the solar r- and s-processes (i.e., rprocess=1 adds the solar r-process amounts
        of each isotope, so to recover the solar isotopic ratios input rprocess=1 and 
        sprocess=1)

    input:
        rprocess : (float) scale factor for the amount of solar pure r-process isotopes
        sprocess : (float) scale factor for the amount of solar pure s-process isotopes

    outupt:
        iso_dict:  (dict) dictionary of isotopes and their relative fractions for a given
            mix of neutron capture processes
    '''


    isotope_data = np.genfromtxt(filepath, dtype=None, names=True, encoding=None, skip_header=2)

    z = isotope_data['z']
    isotopes = isotope_data['isotope']

    r_isotope_abund = isotope_data['nsun'] * isotope_data['rfrac']
    s_isotope_abund = isotope_data['nsun'] * isotope_data['sfrac']

    isotope_abund = (r_isotope_abund * rprocess) + (s_isotope_abund * sprocess)

    iso_dict = calculate_isotope_fractions(z, isotopes, isotope_abund)
    return(iso_dict)


def iso_mix_dict(mix='solar'):
    '''
    Calculates the isotopic ratios of predefined r-/s-process mixes.  Options are to
        calculate the isotope fractions of a pure r-process, pure s-process, or solar mix
        and return these as a dictionary.

    input:
        mix : (str) options are "r", "s", or "solar"

    outupt:
        iso_dict:  (dict) dictionary of isotopes and their relative fractions for a given
            mix of neutron capture processes
    '''

    if mix not in ['solar', 'r', 's']:
        raise ValueError('This isotopic mix is not currently available, '
                         'try using create_iso_dict to specify r- and s-process'
                         'amounts directly')

    if mix == 'r':
        return(create_iso_dict(rprocess=1.0, sprocess=0.0))
    if mix == 's':
        return(create_iso_dict(rprocess=0.0, sprocess=1.0))
    if mix == 'solar':
        return(create_iso_dict(rprocess=1.0, sprocess=1.0))
