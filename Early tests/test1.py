
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.coordinates


# open fits file
file_name = 'data/AS209.nsc.robust_0.5.image.fits'
hdul = fits.open(file_name)




# header
hdr = hdul[0].header
hdr1 = hdul[1].header

print(hdr)
print(list(hdr.keys()))
for card in (hdr.cards):
    print(f"{card.keyword:<8} | {card.value:<20} | {card.comment}")




# Image
# Depending on the data, it may or may not include I, Q, U, V
# data has dimensions of [stokes, frequency, x-axis, y-axis]
data = hdul[0].data
data1 = hdul[1].data
print(data1.shape)
print(data1)

print(data1[0][0])
print(data1[0][1])
print(data1[0][2])
# Beam info
bmaj = hdr1['ttype1']  # [arcsec] -- does this need to be changed? I think it's already in arcsec according to data table
bmin = hdr1['ttype2']  # [arcsec]
bpa = hdr1['ttype3'] # [deg]
print(bmaj)

# Image Coordinates
nx = hdr['naxis2'] # width of x dimension
dx = hdr['cdelt2'] * 3600 # [arcsec] --cdelt refers to scale/step size per pixel
x = (np.arange(nx) - (hdr['crpix2']-1) ) * dx
dec = (np.arange(nx) - (hdr['crpix2']-1)) * hdr['cdelt2']  + hdr['crval2'] # are these mislabled? I thought dec was usually w/ x axis
# then conver to arc sec


ny = hdr['naxis1']
dy = hdr['cdelt1'] * 3600 # [arcsec]
y = (np.arange(ny) - (hdr['crpix1']-1) ) * dy
ra = (np.arange(ny) - (hdr['crpix1']-1)) * hdr['cdelt1'] + hdr['crval1']

#center image


# frequency [Hz]
nf = hdr['naxis3']
f = (np.arange(nf) - (hdr['crpix3'] - 1) ) * hdr['cdelt3'] + hdr['crval3']

# Stokes parameters in units of Jy/beam
I = data[0,0,:,:]
Q = data[1,0,:,:]
U = data[2,0,:,:]
V = data[3,0,:,:]

I_arr = np.array(I)

print(I_arr.shape)
def plot_q():
    ax = plt.gca()
    im = ax.contourf(dec, ra, Q,cmap='viridis', levels=50)
    ax.set_xlabel('x (ra)')
    ax.set_ylabel('y (dec)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Stokes Q [Jy/beam]')

    #plt.xlim(.0004,.001)
    #plt.ylim(-.0095,-.0088)
    plt.show()

def plot_i():
    ax = plt.gca()
    im = ax.contourf(dec, ra, I.T,cmap='viridis', levels=50)
    ax.set_xlabel('x (Dec)')
    ax.set_ylabel('y (RA)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Stokes I [Jy/beam]')

    #plt.xlim(-2,2)
    #plt.ylim(-2,2)
    plt.show()
    plt.close()

def plot_u():
    ax = plt.gca()
    im = ax.contourf(dec, ra, U.T, cmap='viridis', levels=50)
    ax.set_xlabel('x (Dec)')
    ax.set_ylabel('y (RA)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Stokes U [Jy/beam]')

    #plt.xlim(-2, 2)
    #plt.ylim(-2, 2)
    plt.show()



#if __name__ == "__main__":
    # plot_i()
    # plot_q()
    # plot_u()





