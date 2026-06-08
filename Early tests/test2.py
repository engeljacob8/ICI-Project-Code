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
from scipy.interpolate import make_interp_spline
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


P_debiased = np.full_like(P_m, np.nan)

pm_values = np.linspace(0,5*noise_P, 50)
p_values = np.linspace(0,I_arr.max(), 100)
results = np.full_like(pm_values, np.nan)

for i in range(len(pm_values)):
    results[i] = np.argmax(
                    p_values / noise_P ** 2
                    * i0((p_values * pm_values[i]) / (noise_P ** 2))
                    * np.exp(-(pm_values[i] ** 2 + p_values ** 2) / (2 * noise_P ** 2))
                )


interp_obj = make_interp_spline(pm_values, results )
P_debiased = np.where(np.isinf(P_m),0,np.where(P_m < 5 * noise_P,interp_obj(P_m), P_simple))

# linear polarization fraction
error_pf = P_simple/I_arr * np.sqrt((noise_P/P_simple)**2 + (noise_I/I)**2)

pf_debiased = np.full_like(I_arr, np.nan)
pf_simple = np.full_like(I_arr, np.nan)

mask1 = (
        (I_arr > 3* noise_I) &
        (P_simple >3*   noise_P)
)
pf_simple[mask1] = P_simple[mask1]/I_arr[mask1]

mask = (
        (I_arr > 3* noise_I) &
        (P_debiased > 3* noise_P)
    )

#pf = P_simple/I_arr
pf_debiased[mask] = P_debiased[mask]/I_arr[mask]



# #polarization angle
# pa_arr = .5 * np.arctan2(U_arr, Q_arr) # this angle is east of north
# #cartesian
# Ux =  100*pf *-np.sin(pa_arr)
# Uy =  100*pf * np.cos(pa_arr)




if __name__ == "__main__":
    pltf.plot_i(I,ra,dec, noise_I,bmaj= bmaj, bmin= bmin, bpa= bpa, U_arr= U_arr, Q_arr= Q_arr, pf= pf_debiased)

    #pltf.plot_i(I,ra, dec, noise_I)
    # pltf.plot_q(Q, ra, dec, I, noise_I)
    # #
    # pltf.plot_u(U, ra, dec, I, noise_I, noise_U= noise_U)
    # #
    pltf.cont_plot_p(P_debiased, I, noise_I, ra, dec)
    pltf.cont_plot_p(P_simple, I, noise_I, ra, dec)
    # #
    # pltf.im_plot_p(P_simple,extent ,noise_I, I,bmaj= bmaj, bmin= bmin, bpa= bpa)
    # #
    pltf.image_plot_pf(pf_simple, extent)
    pltf.image_plot_pf(pf_debiased, extent)


    hdul.close()



