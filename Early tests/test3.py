import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.special import i0
from scipy.interpolate import make_interp_spline
from Plot_Class import Plot
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

if __name__ == "__main__":
    plot1 = Plot.Plot(ra,dec)
    plot1.set_stokes(I, Q, U)
    plot1.set_noise(noise_I, noise_Q, noise_U)
    plot1.set_beam(bmaj, bmin, bpa)

    ax1 = plot1.plot(I_arr, 'Stokes I', 'Jy/beam', beam=True, contour=True, sig_levels=True)


    obs_angles, vra_azi, vdec_azi = plot1.plot_vect_radius(pf_debiased, ax1, comparison=True)
    ax1.set_xlim(3,-3)
    ax1.set_ylim(-3,3)
    plt.show()
    plot1.compare_angles(obs_angles, vra_azi, vdec_azi)