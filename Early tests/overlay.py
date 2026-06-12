import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.special import i0
from scipy.interpolate import make_interp_spline
from Plot_Class import Plot
# open fits file

file_name = 'data/AS209_continuum.fits'
hdul = fits.open(file_name)




# header
hdr = hdul[0].header


print(hdr)
print(list(hdr.keys()))
for card in (hdr.cards):
    print(f"{card.keyword:<8} | {card.value:<20} | {card.comment}")

# Image
# Depending on the data, it may or may not include I, Q, U, V
# data has dimensions of [stokes, frequency, x-axis, y-axis]
data = hdul[0].data


# Beam info
bmaj = hdr['BMAJ']  # [arcsec] -- does this need to be changed? I think it's already in arcsec according to data table
bmin = hdr['BMIN']  # [arcsec]
bpa = hdr['BPA']  # [deg]

# Image Coordinates
#Q is there a better way to center + change to offset
nx = hdr['naxis2'] # width of x dimension
dx = hdr['cdelt2'] * 3600 # [arcsec] --cdelt refers to scale/step size per pixel
x = (np.arange(nx) - (hdr['crpix2']-1) ) * dx
dec = (np.arange(nx) - (hdr['crpix2']-1)) * hdr['cdelt2'] * 3600  #+ hdr['crval2']
#dec= dec + 14.36927778 # center image/offset

ny = hdr['naxis1']
dy = hdr['cdelt1'] * 3600 # [arcsec]
y = (np.arange(ny) - (hdr['crpix1']-1) ) * dy
ra = (np.arange(ny) - (hdr['crpix1']-1)) * hdr['cdelt1'] * 3600 #+ hdr['crval1']
#ra = ra - 252.31355833 # center image/offset

#map coordinates
extent = [
    ra.max(), ra.min(),
    dec.min(), dec.max(),
]


# center at 16:49:15.254, -14:22:09.4)
# conversion for center
#c = SkyCoord('16h49m15.254s', '-14d22m09.4s', frame = 'icrs')


# frequency [Hz]
nf = hdr['naxis3']
f = (np.arange(nf) - (hdr['crpix3'] - 1) ) * hdr['cdelt3'] + hdr['crval3']