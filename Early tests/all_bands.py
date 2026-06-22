import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.special import i0
from scipy.interpolate import make_interp_spline
from Plot_Class import Plot, Plot_updated


def set_band_7():
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
    I = data_I[0, 0, :, :]
    Q = data_Q[0, 0, :, :]
    U = data_U[0, 0, :, :]

    # Beam info
    bmaj = hdr_I['BMAJ'] * 3600  # [arcsec] --
    bmin = hdr_I['BMIN'] * 3600  # [arcsec]
    bpa = hdr_I['BPA']  # [deg]

    # Image Coordinates
    # Q is there a better way to center + change to offset
    nx = hdr_I['naxis2']  # width of x dimension
    dx = hdr_I['cdelt2'] * 3600  # [arcsec] --cdelt refers to scale/step size per pixel
    x = (np.arange(nx) - (hdr_I['crpix2'] - 1)) * dx
    dec = (np.arange(nx) - (hdr_I['crpix2'] - 1)) * hdr_I['cdelt2'] * 3600  # + hdr_I['crval2']

    ny = hdr_I['naxis1']
    dy = hdr_I['cdelt1'] * 3600  # [arcsec]
    y = (np.arange(ny) - (hdr_I['crpix1'] - 1)) * dy
    ra = (np.arange(ny) - (hdr_I['crpix1'] - 1)) * hdr_I['cdelt1'] * 3600  # + hdr_I['crval1']

    #noise params from carta image
    noise_I = 8.684837817e-5 # Jy/beam
    noise_Q = 2.58228064872e-5
    noise_U = 2.776988582812e-5
    noise_P = np.sqrt(max((noise_Q**2 + noise_U**2)/2,0))

    I_arr = np.array(I)
    U_arr = np.array(U)
    Q_arr = np.array(Q)

    plot1 = Plot_updated.Plot_updated(ra,dec)

    plot1.set_stokes(I, Q, U)
    plot1.set_noise(noise_I, noise_Q, noise_U)
    plot1.set_beam(bmaj, bmin, bpa)
    plot1.set_band('Band 7')

    return plot1

def set_band_6():
    file_name = 'data/member.uid___A001_X136d_X1d5.AS_209_sci.spw25_27_29_31.mfs.IQUV.manual.pbcor (1).fits'
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
    bmin = data1[0][1]  # [arcsec]
    bpa = data1[0][2]  # [deg]

    nx = hdr['naxis2']  #
    dec = (np.arange(nx) - (hdr['crpix2'] - 1)) * hdr['cdelt2'] * 3600  # + hdr['crval2']


    ny = hdr['naxis1']
    ra = (np.arange(ny) - (hdr['crpix1'] - 1)) * hdr['cdelt1'] * 3600  # + hdr['crval1']


    # Stokes parameters in units of Jy/beam
    I = data[0, 0, :, :]
    Q = data[1, 0, :, :]
    U = data[2, 0, :, :]

    # noise params from carta image
    noise_I = 4.137296665e-4  # Jy/beam
    noise_Q = 1.821527e-5
    noise_U = 1.557667505766e-5


    plot1 = Plot_updated.Plot_updated(ra, dec)

    plot1.set_stokes(I, Q, U)
    plot1.set_noise(noise_I, noise_Q, noise_U)
    plot1.set_beam(bmaj, bmin, bpa)
    plot1.set_window(4)
    plot1.set_band('Band 6')

    return plot1

def set_band_3():
    file_name = 'data/AS209.nsc.robust_0.5.image.fits'
    hdul = fits.open(file_name)

    # header
    hdr = hdul[0].header

    # Image
    # Depending on the data, it may or may not include I, Q, U, V
    # data has dimensions of [stokes, frequency, x-axis, y-axis]
    data = hdul[0].data
    data1 = hdul[1].data

    # Beam info
    bmaj = data1[0][0]  # [arcsec] -- does this need to be changed? I think it's already in arcsec according to data table
    bmin = data1[0][1]  # [arcsec]
    bpa = data1[0][2]  # [deg]


    nx = hdr['naxis2']
    dec = (np.arange(nx) - (hdr['crpix2'] - 1)) * hdr['cdelt2'] * 3600  # + hdr['crval2']


    ny = hdr['naxis1']
    ra = (np.arange(ny) - (hdr['crpix1'] - 1)) * hdr['cdelt1'] * 3600  # + hdr['crval1']


    # Stokes parameters in units of Jy/beam
    I = data[0, 0, :, :]
    Q = data[1, 0, :, :]
    U = data[2, 0, :, :]


    # noise params from carta image
    noise_I = 4.939872629984e-5  # Jy/beam
    noise_Q = 8.615981722216e-6
    noise_U = 8.383717527272e-6

    plot1 = Plot_updated.Plot_updated(ra, dec)
    plot1.set_stokes(I, Q, U)
    plot1.set_noise(noise_I, noise_Q, noise_U)
    plot1.set_beam(bmaj, bmin, bpa)
    plot1.set_band('Band 3')

    return plot1



if __name__ == "__main__":
    band7 = set_band_7()
    band6 = set_band_6()
    band3 = set_band_3()
    plot1 = Plot_updated.Plot_updated(band7.ra, band7.dec)
    plot1.set_band('')

    bands = [band7, band6, band3]

    ax = plt.gca()

    plot1.plot_spectrum(bands, ax)


    # ax = band3.function_of_radius(ax,plot_t=True)
    # ax = band7.function_of_radius(ax,plot_t=True)
    # ax = band6.function_of_radius(ax,plot_t=True)
    # plt.show()
    #
    # ax1 = plt.gca()
    #
    # ax1 = band3.function_of_radius(ax1, plot_t=False)
    # ax1 = band7.function_of_radius(ax1, plot_t=False)
    # ax1  = band6.function_of_radius(ax1, plot_t=False)
    # plt.show()





