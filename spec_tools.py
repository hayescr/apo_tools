import numpy as np
import requests
import os
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline



BASE_DIRECTORY = 'https://dr16.sdss.org/sas/dr16/apogee/spectro/aspcap/'
BASE_DIRECTORY_DR17 = 'https://dr17.sdss.org/sas/dr17/apogee/spectro/aspcap/'

library_dict = {
    'synspec_rev1': ['dr17', 'synspec_rev1', 'syn_rev1'],
    'synspec': ['synspec', 'syn'],
    'synspec_lte': ['synspec_lte', 'syn_lte'],
    'turbo20': ['turbo20', 'turbospec', 'turbo', 'ts',
                'turbospec_sph', 'turbo_sph', 'ts_sph'],
    'turbo20_pp': ['turbo20_pp', 'turbospec_pp', 'turbo_pp', 'ts_pp'],
}


class AspcapStar:
    '''
    A class to easily access the aspcapStar 1D, pseudo-normalized spectrum of a
    that has been observed by APOGEE
    '''

    def __init__(self, telescope, field, apogee_id, version='dr16'):
        self.telescope = telescope
        self.field = field
        self.apogee_id = apogee_id
        self.version = 'none'
        if version == 'dr16':
            self.get_dr16()
        else:
            pass

    def get_dr16(self):
        '''
        Gets the DR16 aspcapStar spectrum and fit of a star in APOGEE
        DR16, and saves the wavelength and flux of the spectrum as a Spectrum
        object.
        '''

        self.dr16_obs = Spectrum(
            *self.get_spectrum(BASE_DIRECTORY, 'r12', 'l33', 1))
        self.dr16_fit = Spectrum(
            *self.get_spectrum(BASE_DIRECTORY, 'r12', 'l33', 3))

    def get_spectrum(self, input_directory, reduction, library, extension,
                     path='', save_file=False, **kwargs):
        '''
        Opens an aspcapStar spectrum file from the input directory, reduction
        version and library to the specified path.
        '''

        local_filepath = path + \
            'aspcapStar-{}-{}.fits'.format(reduction, self.apogee_id)

        if not os.path.isfile(local_filepath):
            self.download_spectrum(
                input_directory, reduction, library, local_filepath, **kwargs)

        wave, data = readspec(local_filepath, extension)

        data_m = np.ma.masked_where(data == 0, data)
        data_m = data
        data_m[data_m == 0] = np.nan

        if not save_file:
            os.remove(local_filepath)

        return wave, data_m

    def get_summarytable(self, input_directory, reduction, library, extension,
                         path='', save_file=False, **kwargs):
        '''
        Opens an aspcapStar spectrum file from the input directory, reduction
        version and library to the specified path and returns the summary
        table.
        '''

        local_filepath = path + \
            'aspcapStar-{}-{}.fits'.format(reduction, self.apogee_id)

        if not os.path.isfile(local_filepath):
            self.download_spectrum(
                input_directory, reduction, library, local_filepath, **kwargs)

        with fits.open(local_filepath) as file:
            summary = file[extension].data

        if not save_file:
            os.remove(local_filepath)

        return summary

    def download_dr16(self, path='', **kwargs):
        '''
        Downloads the DR16 aspcapStar spectrum and saves it to the input
        directory (local directory if no path is given).
        '''
        self.download_spectrum(BASE_DIRECTORY, 'r12',
                               'l33', path=path, **kwargs)

    def download_spectrum(self, input_directory, reduction, library,
                          local_filepath, **kwargs):
        '''
        Downloads an aspcapStar spectrum file from the input directory,
        reduction version and library to the specified path.  If
        authentification is required submit requests keywords through **kwargs.
        '''

        fileurl = input_directory + '{}/{}/{}/{}/aspcapStar-{}-{}.fits'.format(
            reduction, library, self.telescope, self.field, reduction,
            self.apogee_id)

        filepath = requests.get(fileurl, **kwargs)
        filepath.raise_for_status()

        with open(local_filepath, 'wb') as file:
            file.write(filepath.content)


class AspcapStarDR17(AspcapStar):
    def __init__(self, telescope, field, apogee_id, version='dr17'):
        super().__init__(telescope, field, apogee_id, version)
        for library, version_options in library_dict.items():
            if version in version_options:
                self.get_version(library)

    def get_version(self, library):
        self.obs = Spectrum(*self.get_spectrum(BASE_DIRECTORY_DR17, 'dr17',
                                               library, 1, save_file=True))
        self.fit = Spectrum(*self.get_spectrum(BASE_DIRECTORY_DR17, 'dr17',
                                               library, 3, save_file=True))
        self.err = Spectrum(*self.get_spectrum(BASE_DIRECTORY_DR17, 'dr17',
                                               library, 2, save_file=True))
        self.info_table = self.get_summarytable(BASE_DIRECTORY_DR17, 'dr17',
                                                library, 4)


class Spectrum:
    def __init__(self, wavelength, flux, variance=None):
        self.wavelength = wavelength
        self.flux = flux
        if variance is not None:
            self.variance = variance

    def vac_to_air(self):
        self.wavelength_air = air_conversion(self.wavelength)

    def vel_shift(self, velocity):
        return velocity_shift(self.wavelength, velocity)

    def interpolate_spectrum(self, new_wavelength):
        ip = interpolate.InterpolatedUnivariateSpline(self.wavelength,
                                                      self.flux,
                                                      k=3, ext='zeros')
        self.wavelength = new_wavelength
        self.flux = ip(new_wavelength)



def readspec(filename, extension):
    '''
    A function designed for easy reading of fits 1D spectra files. Filenames
    can either be input as a local fits file or as a URL to files at a URL
    link.

    Returns the wavelength and flux of the 1D spectrum.
    '''

    with fits.open(filename) as file:
        # file = fits.open(filename)
        flux = file[extension].data
        header = file[extension].header

        ctype = header['CTYPE1']
        cdelt = header['CDELT1']
        crpix = header['CRPIX1']
        cinit = header['CRVAL1']
        naxis1 = header['NAXIS1']

        if ctype in ['LINEAR', 'WAVELENGTH', 'AWAV']:
            wavelength = (np.arange(naxis1) + crpix) * cdelt + cinit

        if ctype == 'LOG-LINEAR':
            wavelength = np.power(10., (np.arange(naxis1)) * cdelt + cinit)

    return wavelength, flux


def vac_spec(filename, extension=1):
    '''
    Opens a vacuum wavelength spectrum and masks pixels with 0 flux.  Returns
    the wavelength and flux as a numpy-like array.
    '''

    wave, data = readspec(filename, extension)

    data_m = np.ma.masked_where(data == 0, data)

    return wave, data_m


def air_spec(filename, extension=1):
    '''
    Opens a vacuum wavelength spectrum, converts its wavelengths to air
    following Shetrone et al. (2015) and masks pixels with 0 flux.  Returns the
    wavelength and flux as a numpy-like array.
    '''

    wave, data = readspec(filename, extension)
    wave_air = air_conversion(wave)
    data_m = np.ma.masked_where(data == 0, data)

    return wave_air, data_m


def air_conversion(wave):
    '''
    Converts a vacuum wavelength spectrum to air following Shetrone et al.
    (2015).
    '''

    a = 0.0
    b1 = 5.792105e-2
    b2 = 1.67917e-3
    c1 = 238.0185
    c2 = 57.362

    wave = wave / 10000.

    air_conv = a + (b1 / (c1 - 1 / (wave**2.))) + \
        (b2 / (c2 - 1 / (wave**2.))) + 1

    wave_air = wave / air_conv
    wave_air = wave_air * 10000.

    return wave_air


def velocity_shift(wavelength, velocity):
    c = 2.99792458e5
    return wavelength * (1 + (velocity / c))
