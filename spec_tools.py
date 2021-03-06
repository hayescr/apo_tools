import numpy as np
import requests
import os
from astropy.io import fits


BASE_DIRECTORY = 'https://dr16.sdss.org/sas/dr16/apogee/spectro/aspcap/'


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
                     save_file=False, **kwargs):
        '''
        Opens an aspcapStar spectrum file from the input directory, reduction
        version and library to the specified path.
        '''

        local_filepath = self.download_spectrum(
            input_directory, reduction, library, **kwargs)

        wave, data = readspec(local_filepath, extension)

        data_m = np.ma.masked_where(data == 0, data)
        data_m = data
        data_m[data_m == 0] = np.nan

        if not save_file:
            os.remove(local_filepath)

        return wave, data_m

    def download_dr16(self, path='', **kwargs):
        '''
        Downloads the DR16 aspcapStar spectrum and saves it to the input
        directory (local directory if no path is given).
        '''
        self.download_spectrum(BASE_DIRECTORY, 'r12',
                               'l33', path=path, **kwargs)

    def download_spectrum(self, input_directory, reduction, library, path='', **kwargs):
        '''
        Downloads an aspcapStar spectrum file from the input directory, reduction
        version and library to the specified path.  If authentification is
        required submit requests keywords through **kwargs.
        '''

        fileurl = input_directory + '{}/{}/{}/{}/aspcapStar-{}-{}.fits'.format(
            reduction, library, self.telescope, self.field, reduction,
            self.apogee_id)

        filepath = requests.get(fileurl, **kwargs)
        filepath.raise_for_status()

        local_filepath = path + \
            'aspcapStar-{}-{}.fits'.format(reduction, self.apogee_id)

        with open(local_filepath, 'wb') as file:
            file.write(filepath.content)

        return local_filepath


class Spectrum:
    def __init__(self, wavelength, flux):
        self.wavelength = wavelength
        self.flux = flux

    def vac_to_air(self):
        self.wavelength_air = air_conversion(self.wavelength)


def readspec(filename, extension):
    '''
    A function designed for easy reading of fits 1D spectra files. Filenames can
    either be input as a local fits file or as a URL to files at a URL link.

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

        if ctype == 'LINEAR' or ctype == 'WAVELENGTH':
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
    (2015) and masks pixels with 0 flux.  Returns the wavelength and flux as a
    numpy-like array.
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
