from turtledemo.chaos import plot

from matplotlib import pyplot as plt
import numpy as np
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

    def set_beam(self, bmaj, bmin, bpa):
        self.bmaj = bmaj
        self.bmin = bmin
        self.bpa = bpa

    def set_noise(self, noise_I, noise_Q,noise_U):
        self.noise_I = noise_I
        self.noise_Q = noise_Q
        self.noise_U = noise_U
        self.noise_P = noise_P = np.sqrt(max((noise_Q**2 + noise_U**2)/2,0))

    def set_stokes(self, I, Q, U):
        self.I = np.array(I)
        self.Q = np.array(Q)
        self.U = np.array(U)



    def plot_beam(self,ax):
        my_beam = Beam(self.bmaj * u.arcsec, self.bmin * u.arcsec, self.bpa * u.degree)
        y_cen_pix, x_cen_pix = -2, 2
        pixscale = u.arcsec
        ellipse_artist = my_beam.ellipse_to_plot(x_cen_pix, y_cen_pix, pixscale, facecolor='black')
        ax.add_artist(ellipse_artist)


    def plot (self, value, value_name, unit, beam=False,contour=False, sig_levels=False, ax=None, color='viridis'):
        if ax is None:
            ax = plt.gca()
        if contour:
            im = ax.contourf(self.ra, self.dec, value, cmap=color, levels=50)
        else:
            #using convention of east to left, north to up
            extent = [
                self.ra.max(), self.ra.min(),
                self.dec.min(), self.dec.max(),
            ]
            im = ax.imshow( value,cmap=color, extent= extent)

        ax.set_xlabel(r'$\Delta$RA (arcsec)')
        ax.set_ylabel(r'$\Delta$Dec (arcsec)')
        ax.set_aspect('equal')

        if beam:
            self.plot_beam(ax)

        if sig_levels:
            contours = ax.contour(self.ra, self.dec, self.I, colors='white', linewidths=.6,
                                  levels=[3 * self.noise_I, 10 * self.noise_I, 25 * self.noise_I, 50 * self.noise_I, 100 * self.noise_I,
                                          200 * self.noise_I,
                                          325 * self.noise_I, 500 * self.noise_I, 1000 * self.noise_I])
            ax.clabel(contours, inline=1, fontsize=10)

        cbar = plt.colorbar(im)
        cbar.set_label(f'{value_name} [{unit}]')

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
        # .3 arcsec for 1% polarization
        d_ra = .3 * pf_sample * 100 * np.sin(pa_arr) #dec on x ax
        d_dec = .3 * pf_sample * 100 * np.cos(pa_arr)

        comparison_angle = np.atan2(vra_azi,vdec_azi)
        vra = np.sin(comparison_angle)
        vdec = np.cos(comparison_angle)
        #ax.quiver(ra_sample, dec_sample, -vy_sky, vx_sky, pivot='middle',color='white', scale=5, scale_units='xy',headwidth=1e-10,
                  #headlength=1e-10, headaxislength=1e-10, width=0.005)

        q = ax.quiver(ra_sample, dec_sample, d_ra, d_dec, pivot='middle', angles='xy', color='red', scale=1, scale_units='xy', headwidth=1e-10,
                  headlength=1e-10, headaxislength=1e-10, width=0.005)
        ax.quiverkey(q, X=0.85, Y=1.05, U=.3, label='1 % Polarization', labelpos='E')
        if comparison:
            ax.quiver(ra_sample, dec_sample, vra, vdec, pivot='middle',angles='xy',color='white', scale=5, scale_units='xy',headwidth=1e-10,
                 headlength=1e-10, headaxislength=1e-10, width=0.005)
        return pa_arr, vra_azi, vdec_azi
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
        spacing = 10 #simple start -
        phi = np.linspace(0, 2*np.pi, spacing) # for now, 10
        #radius of 1 arc sec - change later for multiple radii at different points with correct 1/2 beam spacing
        #self.bmaj, self.bmaj + .5
        radii = np.array([self.bmaj/2 ])
        ra_disk = (radii[:,None] * np.cos(phi)).ravel()
        dec_disk = (radii[:,None] * np.sin(phi)).ravel()

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


    def compare_angles(self, obs_angle, vra_azi, vdec_azi):
        exp_angle = np.atan2(vra_azi, vdec_azi)
        diff_angle = np.degrees(obs_angle - exp_angle)
        ax = plt.gca()
        ax.hist(diff_angle)
        ax.set_xlabel('Difference of Observed and Expected Angle')
        ax.set_ylabel('Frequency')
        plt.show()

