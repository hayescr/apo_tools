# apo_tools

This is a repository of various tools that I use to work with, analyze and visualize APOGEE data.

## aspcap_star_example_notebook.ipynb

aspcap_star_example_notebook.ipynb is a Jupyter notebook that gives a walkthrough for using the AspcapStar class in spec_tools.py, and the window plotting functions in spec_plotting.py to plot download and plot APOGEE spectra in a variety of different ways.

## galcoords.py

galcoords.py contains functions to convert between observed values in equatorial coordinates, such as RA, Dec, distance, radial velocity and proper motion to galactic coordinates, or to convert galactic coordinates to the Sagittarius coordinate system as defined by Majewski et al. 2003 and Law and Majewski 2010.

Methods of interest:
1. calculate_gal_vel - calculates Galactocentric velocities and positions
2. calculate_sgr_system - calculates Galactocentric velocities and positions in the Sgr coordinate system
3. calculate_sgr_lb - calculates Sgr sky coordinates lambda and beta


Usage with distances in parsecs (set up for Gaia proper motions which are reported in mas/yr):
```python
from galcoords import Galcoords
new_coords = Galcoords(ra=ra, dec=dec, l=glon,
                       b=glat, pm_ra=gaia_pmra / 1000.,
                       pm_ra_e=gaia_pmra_error / 1000.,
                       pm_dec=gaia_pmdec / 1000.,
                       pm_dec_e=gaia_pmdec_error / 1000.,
                       v_rad=v_rad, v_rad_e=v_rad_error,
                       dist=dist,
                       sigma_dist=dist_error, use_dist=True)
new_coords.calculate_gal_vel()
```


Usage with parallaxes (set up for Gaia parallaxes and proper motions which are reported in mas and mas/yr):

```python
from galcoords import Galcoords
new_coords = Galcoords(ra=ra, dec=dec, l=glon,
                       b=glat, pm_ra=gaia_pmra / 1000.,
                       pm_ra_e=gaia_pmra_error / 1000.,
                       pm_dec=gaia_pmdec / 1000.,
                       pm_dec_e=gaia_pmdec_error / 1000.,
                       v_rad=v_rad, v_rad_e=v_rad_error,
                       parallax=gaia_parallax / 1000.,
                       parallax_e=gaia_parallax_error/ 1000.)
new_coords.calculate_gal_vel()
```

Once these methods have been used you can access the positions and velocities as object attributes (new_coord.xxxxx) and their estimated uncertainties (new_coord.sigma_xxxxx):

- dist - distance (kpc)
- U, V, W - heliocentric cartesian velocities (km/s)
- x, y, z, vx, vy, vz - Galactocentric cartesian positions (kpc) and velocities (km/s)
- R_cyn, phi, vR, vphi - Galactocentric cylindrical positions and velocities (km/s)

Calculating Sgr system coordinates, positions and velocities give similar values but adding an "s" to the end (except r_cys for the cylindrical position), and lambda_sun and beta_sun the solar centered Sgr stream coordinates and lambda_gc and beta_gc, the Galactocentric Sgr stream coordinates.

## spec_tools.py

spec_tools.py contains functions to work with APOGEE spectra, such as opening or downloading spectra from the web and tools to read APOGEE spectra files, mask empty pixels, and convert between vacuum and air wavelengths.

## spec_plotting.py

spec_plotting.py contains functions useful for plotting APOGEE spectra, such as functions to plot the elemental abundance windows in the windows folder.

## cno_combine.py

cno_combine.py provides a function to calculate the combined abundance of C and N - [(C+N)/Fe], and C, N, and O - [(C+N+O)/Fe].  Useful for giants where C and N are not conserved and altered due to dredge up, but their sum is nearly conserved and C, N, and O should be conserved together.

## abund_utils.py

abund_utils.py provides some basic utiilities for working with abundances, such as dictionaries of solar abundances, and functions for converting between atomic symbol and number
