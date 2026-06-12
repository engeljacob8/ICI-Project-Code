from turtledemo.chaos import plot

from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from radio_beam import Beam
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator as rgi

class Plot:
    """Class for plotting polarization plots with RA as x-axis and DEC as y-axis"""

    def __init__(self, ra, dec):
        self.I = None
        self.Q = None
        self.U = None
        self.bmaj = None
        self.bmin = None
        self.bpa = None
        self.noise_I = None
        self.noise_Q = None
        self.noise_U = None
        self.noise_P = None
        self.angle_incl = 34.88
        self.angle_pa = 85.8
        self.ra = ra
        self.dec = dec
        self.window = 3

    def set_beam(self, bmaj, bmin, bpa):
        self.bmaj = bmaj
        self.bmin = bmin
        self.bpa = bpa

    def set_band(self,band):
        self.band = band

    def set_window(self, window):
        self.window = window

    def set_noise(self, noise_I, noise_Q,noise_U):
        self.noise_I = noise_I
        self.noise_Q = noise_Q
        self.noise_U = noise_U
        self.noise_P = noise_P = np.sqrt(max((noise_Q**2 + noise_U**2)/2,0))

    def set_stokes(self, I, Q, U):
        self.I = np.array(I)
        self.Q = np.array(Q)
        self.U = np.array(U)



    def plot_beam(self,ax, bmaj=None, bmin=None, bpa=None):
        if bmaj is None:
            my_beam = Beam(self.bmaj * u.arcsec, self.bmin * u.arcsec, -self.bpa * u.degree)
            y_cen_pix, x_cen_pix = -self.window + 1, self.window - 1
            color= 'black'
        else:
            my_beam = Beam(bmaj * u.arcsec, bmin * u.arcsec, -bpa * u.degree)
            y_cen_pix, x_cen_pix = -self.window + 1, -self.window + 1
            color = 'white'
            #ax.annotate(text='100 AU', xy=(1, 1.8), xytext=(-self.window + 1, -self.window + 1), textcoords='data', fontsize=7,
                      #  color='black', weight='bold')

        pixscale = u.arcsec
        ellipse_artist = my_beam.ellipse_to_plot(x_cen_pix, y_cen_pix, pixscale, facecolor=color)
        ax.add_artist(ellipse_artist)


    def plot (self, value, value_name, unit,overlay_extent=None, overlay=False, beam=False,contour=False, sig_levels=False, ax=None,au=False ,color='viridis'):
        if ax is None:
            ax = plt.gca()
        if contour:
            if unit == '%':
                levels = [0,.3,.6,.9,1.2,1.5]
                im = ax.contourf(self.ra, self.dec, value * 100, cmap=color, levels=levels, extend = 'max')
            elif value_name == 'Stokes Q':
                contours1 = ax.contour(self.ra, self.dec, value, colors='red', linewidths=.5, levels=[3 * self.noise_Q ])
                contours2 = ax.contour(self.ra, self.dec, value, colors='blue', linewidths=.5, levels=[-3 * self.noise_Q])
                im = ax.contourf(self.ra, self.dec, value, cmap=color, levels=50)
            elif value_name == 'Stokes U':
                contours1 = ax.contour(self.ra, self.dec, value, colors='red', linewidths=.5, levels=[3 * self.noise_U ])
                contours2 = ax.contour(self.ra, self.dec, value, colors='blue', linewidths=.5, levels=[-3 * self.noise_U])
                im = ax.contourf(self.ra, self.dec, value, cmap=color, levels=50)
            else:
                im = ax.contourf(self.ra, self.dec, value, cmap=color, levels=50)
        elif overlay:
            im = ax.imshow(value, cmap=color, extent= overlay_extent)
        else:
            #using convention of east to left, north to up
            extent = [
                self.ra.max(), self.ra.min(),
                self.dec.min(), self.dec.max(),
            ]
            im = ax.imshow(value,cmap=color, extent= extent)

        ax.set_xlabel(r'$\Delta$RA (arcsec)')
        ax.set_ylabel(r'$\Delta$Dec (arcsec)')
        ax.set_aspect('equal')

        if overlay:
            ax.set_title(f'Polarization at band {self.band} + 1.2mm Dust')
        else:
            ax.set_title(f'{value_name} {self.band}')

        if beam:
            self.plot_beam(ax)
        if au:
            ax.quiver(1,-self.window+1, 0.82644628099,0,  pivot='middle', angles='xy', color='black', scale=1, scale_units='xy', headwidth=1e-10,
                  headlength=1e-10, headaxislength=1e-10, width=0.005)
            ax.annotate(text='100 AU', xy=(1,1.8), xytext=(1.4,-self.window+1.1), textcoords='data', fontsize=7, color='black', weight='bold')

        if sig_levels:
            contours = ax.contour(self.ra, self.dec, self.I, colors='white', linewidths=.6, alpha=.5,
                                  levels=[3 * self.noise_I, 10 * self.noise_I, 25 * self.noise_I, 50 * self.noise_I, 100 * self.noise_I,
                                          200 * self.noise_I,
                                          325 * self.noise_I, 500 * self.noise_I, 1000 * self.noise_I])


        cbar = plt.colorbar(im)
        cbar.set_label(f'{value_name} [{unit}]')

        ax.set_xlim(self.window, -self.window)
        ax.set_ylim(-self.window, self.window)
        #change later - but base limit and axis switch
        # plt.xlim(3, -3)
        # plt.ylim(-3, 3)
        return ax



    def plot_vect_radius(self,pf,ax, comparison=False):
        #sample points along radius
        ra_sample, dec_sample, vra_azi, vdec_azi = self.create_radius()

        # instantiate rectangular grid interpolator
        Q_interp = rgi((self.dec,self.ra), self.Q, bounds_error=False, fill_value=np.nan) #(dec, ra)
        U_interp = rgi((self.dec,self.ra), self.U, bounds_error=False, fill_value=np.nan)
        pf_interp = rgi((self.dec,self.ra), pf, bounds_error=False, fill_value=np.nan)

        points = np.column_stack((dec_sample, ra_sample))
        # interpolate from given points
        Q_sample = Q_interp(points)
        U_sample = U_interp(points)
        pf_sample = pf_interp(points)

        Q_sample = Q_sample.reshape(dec_sample.shape)
        U_sample = U_sample.reshape(dec_sample.shape)
        pf_sample = pf_sample.reshape(dec_sample.shape)

        # calculate polarization angle
        pa_arr = .5 * np.arctan2(U_sample, Q_sample)  # this angle is east of north
        # cartesian
        # arcsec for 1% polarization
        scale = .5 #arcsec/percent
        d_ra = scale * pf_sample * 100 * np.sin(pa_arr) #dec on x ax
        d_dec = scale * pf_sample * 100 * np.cos(pa_arr)
        pa_arr = np.arctan2(d_ra, d_dec)

        comparison_angle = np.arctan2(vra_azi,vdec_azi)
        vra = np.sin(comparison_angle)
        vdec = np.cos(comparison_angle)


        q = ax.quiver(ra_sample, dec_sample, d_ra, d_dec, pivot='middle', angles='xy', color='red', scale=1, scale_units='xy', headwidth=1e-10,
                  headlength=1e-10, headaxislength=1e-10, width=0.005)
        ax.quiverkey(q, X=0.85, Y=1.05, U=scale, label=f'{scale} % Polarization', labelpos='E')
        if comparison:
            ax.quiver(ra_sample, dec_sample, vra, vdec, pivot='middle',angles='xy',color='white', scale=5, scale_units='xy',headwidth=1e-10,
                 headlength=1e-10, headaxislength=1e-10, width=0.005)
            ax.set_title('')


        return pa_arr, comparison_angle
        #debug check
        # ax.quiver(
        #     [0],
        #     [0],
        #     [1],
        #     [0],
        #     color='yellow',
        #     angles='xy',
        #     scale=1,
        #     scale_units='xy'
        # )
        # ax.quiver(
        #     [0],
        #     [0],
        #     [0],
        #     [1],
        #     color='cyan',
        #     angles='xy',
        #     scale=1,
        #     scale_units='xy'
        # )



    def create_radius(self ):

        radii = np.array([self.bmaj / 4, .75 * self.bmaj, 1.25 * self.bmaj, 1.75 * self.bmaj])
        arc_spacing = self.bmaj / 2

        ra_list = []
        dec_list = []

        for r in radii:
            #approximation for circular distances of 1/2 bmaj - can be improved with elliptical circum
            spacing = max(4, int(np.ceil(2 * np.pi * r / arc_spacing)))
            phi = np.linspace(0, 2 * np.pi, spacing, endpoint=False)

            ra_list.append(r * np.cos(phi))
            dec_list.append(r * np.sin(phi))

        ra_disk = np.concatenate(ra_list)
        dec_disk = np.concatenate(dec_list)



        #create azimuthal comparison
        ra_disk_azi = -dec_disk
        dec_disk_azi = ra_disk

        #convert to image coordinates
        ra_im0 = ra_disk
        dec_im0 = dec_disk * np.cos(np.radians(self.angle_incl))

        vra_im =  ra_disk_azi
        vdec_im = dec_disk_azi * np.cos(np.radians(self.angle_incl))

        # rotate - intentionally changed because east to north is swapped
        angle_pa = np.radians(self.angle_pa)
        ra_image = ra_im0 * np.sin(angle_pa) - np.cos(angle_pa) * dec_im0
        dec_image = dec_im0 * np.sin(angle_pa) + np.cos(angle_pa) * ra_im0

        vra_image= vra_im * np.sin(angle_pa) - np.cos(angle_pa) * vdec_im
        vdec_image = vdec_im * np.sin(angle_pa) + np.cos(angle_pa) * vra_im
        return ra_image, dec_image, vra_image, vdec_image


    def compare_angles(self, obs_angle, exp_angle):
        diff_angle = np.degrees(obs_angle - exp_angle)
        diff_angle = (diff_angle + 90) % 180 - 90
        ax = plt.gca()
        ax.hist(diff_angle, alpha = .7, color = 'blue')
        ax.set_xlabel('Difference of Observed and Expected Angle')
        ax.set_ylabel('Frequency')

        mask = np.isfinite(diff_angle)
        data_mean = np.mean(diff_angle[mask])

        ax.axvline(float(data_mean), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {data_mean:.2f}')

        plt.show()


    def plot_overlay(self,file_name):

        hdul = fits.open(file_name)

        # header
        hdr = hdul[0].header
        data = hdul[0].data

        # Beam info
        bmaj = hdr['BMAJ'] *3600 # [arcsec] --
        bmin = hdr['BMIN']  *3600# [arcsec]
        bpa = hdr['BPA']  # [deg]

        # Image Coordinates
        nx = hdr['naxis2']  # width of x dimension
        dx = hdr['cdelt2'] * 3600  # [arcsec] --cdelt refers to scale/step size per pixel
        x = (np.arange(nx) - (hdr['crpix2'] - 1)) * dx
        dec = (np.arange(nx) - (hdr['crpix2'] - 1)) * hdr['cdelt2'] * 3600  # + hdr['crval2']
        # dec= dec + 14.36927778 # center image/offset

        ny = hdr['naxis1']
        dy = hdr['cdelt1'] * 3600  # [arcsec]
        y = (np.arange(ny) - (hdr['crpix1'] - 1)) * dy
        ra = (np.arange(ny) - (hdr['crpix1'] - 1)) * hdr['cdelt1'] * 3600  # + hdr['crval1']
        # ra = ra - 252.31355833 # center image/offset

        # map coordinates
        extent = [
            ra.max(), ra.min(),
            dec.min(), dec.max(),
        ]

        I = data[0, 0, :, :]

        ax = self.plot(I, 'polarization + 1.2mm dust cont', 'Jy/beam', overlay=True, overlay_extent=extent)
        self.plot_beam(ax, bmaj, bmin, bpa)
        return ax







