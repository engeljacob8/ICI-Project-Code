import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.coordinates
from matplotlib.pyplot import subplots
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.special import i0
from scipy.optimize import minimize
import Cleaned_Plot_Methods.Plot_Functions as pltf

# open fits file
file_name = 'data/AS209.nsc.robust_0.5.image.fits'
hdul = fits.open(file_name)




# header
hdr = hdul[0].header
hdr1 = hdul[1].header


# Image
# Depending on the data, it may or may not include I, Q, U, V
# data has dimensions of [stokes, frequency, x-axis, y-axis]
data = hdul[0].data
data1 = hdul[1].data

# Beam info
bmaj = data1[0][0]  # [arcsec] -- does this need to be changed? I think it's already in arcsec according to data table
bmin = data1[0][1]   # [arcsec]
bpa = data1[0][2]  # [deg]

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

# Stokes parameters in units of Jy/beam
I = data[0,0,:,:]
Q = data[1,0,:,:]
U = data[2,0,:,:]
V = data[3,0,:,:]



#noise params from carta image
noise_I = 4.939872629984e-5 # Jy/beam
noise_Q = 8.615981722216e-6
noise_U = 8.383717527272e-6
noise_P = np.sqrt(max((noise_Q**2 + noise_U**2)/2,0))

I_arr = np.array(I)
U_arr = np.array(U)
Q_arr = np.array(Q)

#polarized intensity
P_m = np.sqrt(Q_arr**2 + U_arr**2)

P_simple = np.sqrt(np.maximum(Q_arr**2 + U_arr**2 - noise_P**2,0))

#debiased intensity
#P = sp.symbols('P')

#pdf = P/noise_P**2 * j0((P * P_m)/noise_P**2) * np.exp(- (P_m**2 + P**2)/2* noise_P**2)
#cdf = 1/(np.pi * noise_P**2) * np.exp(-(P**2 + P_m**2)/(2 * noise_P**2)) ((1-P**2/noise_P**2) * sp.integrate(np.exp(((P*P_m)/noise_P**2)*np.cos(theta)), 0, np.pi ) + P*P_m/noise_P**2 * sp.integrate(np.cos(theta)*np.exp((P*P_m*np.cos(theta)/noise_P**2)), 0, np.pi))

P_debiased = np.full_like(P_m, np.nan)

# for iy in range(P_m.shape[0]):
#     for ix in range(P_m.shape[1]):
#         pm = P_m[iy, ix]
#         if not np.isfinite(pm):
#             continue
#
#         if pm <= 5 * noise_P:
#
#             def neg_pdf(P):
#                 P = x[0]
#
#                 pdf = (
#                     P / noise_P ** 2
#                     * i0((P * pm) / (noise_P ** 2))
#                     * np.exp(-(pm ** 2 + P ** 2) / (2 * noise_P ** 2))
#                 )
#                 return -pdf
#
#             res = minimize(neg_pdf, [max(pm, 1e-10)], bounds=[(0, None)])
#
#             P_debiased[iy, ix] = res.x[0]
#         else:
#             P_debiased[iy, ix] = np.sqrt(max(pm**2 - noise_P**2,0))


# linear polarization fraction
noise_pf = P_simple/I_arr * np.sqrt((noise_P/P_simple)**2 + (noise_I/I)**2)
pf = np.full_like(I_arr, np.nan)

mask = (

        (I_arr > 2* noise_I) &
        (P_simple >2*   noise_P)
    )

#pf = P_simple/I_arr
pf[mask] = P_simple[mask]/I_arr[mask]



#polarization angle
pa_arr = .5 * np.arctan2(U_arr, Q_arr) # this angle is east of north
#cartesian
Ux =  100*pf *-np.sin(pa_arr)
Uy =  100*pf * np.cos(pa_arr)




if __name__ == "__main__":
    pltf.plot_i(I,ra,dec, noise_I, bmaj= bmaj, bmin= bmin, bpa= bpa, ux=Ux, uy=Uy)

    #pltf.plot_q(Q, ra, dec, I, noise_I)
    #
    #pltf.plot_u(U, ra, dec, I, noise_I, noise_U= noise_U)
    #
    # pltf.cont_plot_p(P_simple, I, noise_I, extent, ra, dec)
    #
    #pltf.im_plot_p(P_simple,extent ,noise_I, I,bmaj= bmaj, bmin= bmin, bpa= bpa)
    #
    pltf.image_plot_pf(pf, extent)


    hdul.close()



