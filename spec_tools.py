import numpy as np
from astropy.io import fits


def readspec(filename, extension):
    '''
    A function designed for easy reading of fits 1D spectra files. Filenames can
    either be input as a local fits file or as a URL to files at a URL link.

    Returns the wavelength and flux of the 1D spectrum.
    '''

    file = fits.open(filename)
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


def open_aspcap_spec(telescope, field, apogeeid, extension=1):
    '''
    Accesses the aspcapStar 1D, pseudo-normalized spectrum of a star in APOGEE
    DR16, opens the file and returns the wavelength and flux of the spectrum
    (pixels where there is no flux are masked).
    '''

    fileurl = f'https://dr16.sdss.org/sas/dr16/apogee/spectro/aspcap/r12/l33/{telescope}/{field}/aspcapStar-r12-{apogeeid}.fits'
    wave, data = readspec(fileurl, extension)
    data_m = np.ma.masked_where(data == 0, data)

    return wave, data_m


def download_aspcap_spec(telescope, field, apogeeid, path=''):
    '''
    Downloads the aspcapStar 1D pseudo-normalized spectrum of a star in APOGEE
    DR16 and saves it to the input directory (local directory if no path is given)
    if the file is not already present.
    '''

    if not os.path.isfile(path + 'aspcapStar-r12-{}.fits'.format(apogeeid)):
        fileurl = f'https://dr16.sdss.org/sas/dr16/apogee/spectro/aspcap/r12/l33/{telescope}/{field}/aspcapStar-r12-{apogeeid}.fits'
        filepath = requests.get(fileurl)
        open(path + 'aspcapStar-r12-{}.fits'.format(apogeeid),
             'wb').write(filepath.content)


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
    data_m = np.ma.masked_where(data == 0, data)

    return wave_air, data_m
