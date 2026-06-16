import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.special import i0
from scipy.interpolate import make_interp_spline
from Plot_Class import Plot
# open fits file

file_name_I = 'data/member.uid___A001_X1273_Xb73.as209_sci.spw13_15_17_19.cont.I.manual.clean.pbcor.fits'
file_name_Q = 'data/member.uid___A001_X1273_Xb73.as209_sci.spw13_15_17_19.cont.Q.manual.clean.pbcor.fits'
file_name_U = 'data/member.uid___A001_X1273_Xb73.as209_sci.spw13_15_17_19.cont.U.manual.clean.pbcor.fits'
hdul_I = fits.open(file_name_I)
hdul_Q = fits.open(file_name_Q)
hdul_U = fits.open(file_name_U)

file_name_cont = 'data/AS209_continuum.fits'




# header
hdr_I = hdul_I[0].header
hdr_Q = hdul_Q[0].header
hdr_U = hdul_U[0].header




# Image
# Depending on the data, it may or may not include I, Q, U, V
# data has dimensions of [stokes, frequency, x-axis, y-axis]
data_I = hdul_I[0].data
data_Q = hdul_Q[0].data
data_U = hdul_U[0].data

# Stokes parameters in units of Jy/beam
I = data_I[0,0,:,:]
Q = data_Q[0,0,:,:]
U = data_U[0,0,:,:]



# Beam info
bmaj = hdr_I['BMAJ'] * 3600 #[arcsec] --
bmin = hdr_I['BMIN'] * 3600  # [arcsec]
bpa = hdr_I['BPA'] # [deg]

# Image Coordinates
#Q is there a better way to center + change to offset
nx = hdr_I['naxis2'] # width of x dimension
dx = hdr_I['cdelt2'] * 3600 # [arcsec] --cdelt refers to scale/step size per pixel
x = (np.arange(nx) - (hdr_I['crpix2']-1) ) * dx
dec = (np.arange(nx) - (hdr_I['crpix2']-1)) * hdr_I['cdelt2'] * 3600  #+ hdr_I['crval2']

ny = hdr_I['naxis1']
dy = hdr_I['cdelt1'] * 3600 # [arcsec]
y = (np.arange(ny) - (hdr_I['crpix1']-1) ) * dy
ra = (np.arange(ny) - (hdr_I['crpix1']-1)) * hdr_I['cdelt1'] * 3600 #+ hdr_I['crval1']

#map coordinates
extent = [
    ra.max(), ra.min(),
    dec.min(), dec.max(),
]


# center at 16:49:15.254, -14:22:09.4)
# conversion for center
#c = SkyCoord('16h49m15.254s', '-14d22m09.4s', frame = 'icrs')


# frequency [Hz]
nf = hdr_I['naxis3']
f = (np.arange(nf) - (hdr_I['crpix3'] - 1) ) * hdr_I['cdelt3'] + hdr_I['crval3']




#noise params from carta image
noise_I = 8.684837817e-5 # Jy/beam
noise_Q = 2.58228064872e-5
noise_U = 2.776988582812e-5
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
P_debiased = np.where(P_m < 5 * noise_P,interp_obj(P_m), P_simple)

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
pf_debiased = np.where(mask, P_debiased/I_arr, np.nan)


if __name__ == "__main__":
    plot1 = Plot.Plot(ra,dec)

    plot1.set_stokes(I, Q, U)
    plot1.set_noise(noise_I, noise_Q, noise_U)
    plot1.set_beam(bmaj, bmin, bpa)
    plot1.set_band('Band 7')


    ax = plot1.plot_principle_frame('Stokes Q', 'Jy/beam')
    plt.show()

    axI = plot1.plot_principle_frame('Stokes I', 'Jy/beam')
    obs_angles, comp_angles, chi_error = plot1.plot_vect_radius( axI, comparison=True, principle_frame=True)
    plt.show()
    plot1.compare_angles(obs_angles, comp_angles, chi_error)

    plt.show()

    ax1 = plot1.plot(I_arr, 'Stokes I', 'Jy/beam', beam=True, contour=True, sig_levels=True, au=True,axis=True)
    obs_angles_sky, comp_angles1,chi_error1 = plot1.plot_vect_radius( ax1, comparison=True)
    plt.show()

    #orginal comparison

    plot1.compare_angles(obs_angles_sky, comp_angles1,chi_error1)

    plt.show()
    #test
    plot1.test1(obs_angles, obs_angles_sky)
    plot1.test2()
    plot1.test3()
    plot1.run_all_tests()


    # ax2 = plot1.plot(Q_arr, 'Stokes Q', 'Jy/beam', beam=True, contour=True, sig_levels=True, au=True)
    # plt.show()
    # ax3 = plot1.plot(U_arr, 'Stokes U', 'Jy/beam', beam=True, contour=True, sig_levels=True, au=True)
    # plt.show()
    # ax4 = plot1.plot(P_debiased, 'Stokes P', 'Jy/beam', beam=True, contour=True, sig_levels=True, au=True)
    # plt.show()
    # ax5 = plot1.plot(pf_debiased, 'Polarization Fraction', '%', beam=True, contour=True, sig_levels=True, au=True)
    # plt.show()
    # ax6 = plot1.plot(pf_simple, 'Polarization Fraction', '%', beam=True, contour=True, sig_levels=False)
    # plt.show()
    #
    # ax7 = plot1.plot_overlay(file_name_cont)
    # plot1.plot_vect_radius(pf_debiased, ax7, comparison=False)
    # plt.show()
