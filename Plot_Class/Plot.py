from turtledemo.chaos import plot

from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from radio_beam import Beam
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.special import i0
from scipy.interpolate import make_interp_spline

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
        self.set_principle_frame(None,None,None, None)
        self.pf = None

    def set_pf(self,pf):
        self.pf = pf

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

    def set_principle_frame(self,  Q_princ, U_princ, x_princ, y_princ):

        self.Q_princ = Q_princ
        self.U_princ = U_princ
        self.x_princ = x_princ
        self.y_princ = y_princ




    def plot_beam(self,ax, bmaj=None, bmin=None, bpa=None, overlay = False):
        if bmaj is None:
            my_beam = Beam(self.bmaj * u.arcsec, self.bmin * u.arcsec, -self.bpa * u.degree)
            y_cen_pix, x_cen_pix = -self.window + 1, self.window - 1
            color= 'black'
        elif overlay:
            my_beam = Beam(bmaj * u.arcsec, bmin * u.arcsec, -bpa * u.degree)
            y_cen_pix, x_cen_pix = -self.window + 1, -self.window + 1
            color = 'white'
            #ax.annotate(text='100 AU', xy=(1, 1.8), xytext=(-self.window + 1, -self.window + 1), textcoords='data', fontsize=7,
                      #  color='black', weight='bold')
        else:
            my_beam = Beam(bmaj * u.arcsec, bmin * u.arcsec, -bpa * u.degree)
            y_cen_pix, x_cen_pix = -self.window + 1, self.window - 1
            color = 'black'

        pixscale = u.arcsec
        ellipse_artist = my_beam.ellipse_to_plot(x_cen_pix, y_cen_pix, pixscale, facecolor=color)
        ax.add_artist(ellipse_artist)


    def plot (self, value, value_name, unit,overlay_extent=None,overlay=False,
              beam=False,contour=False, sig_levels=False, ax=None,au=False ,color='viridis', axis=False, principle=False):
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



        if axis:
            slope = np.tan(np.radians(90-self.angle_pa))
            ax.axline((0,0),slope=slope, linestyle='--' ,linewidth=0.5, color='black',alpha=.5)
            ax.axline((0, 0), slope= -1/slope, linestyle='--',linewidth=0.5, color='black',alpha=.5)
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


    def plot_principle_frame(self, value_name, unit, ax=None, color='viridis' ):
        if ax is None:
            ax = plt.gca()
        dec_principle, ra_principle, q_principle, u_principle, maj, x = self.principle_frame()

        if value_name == 'Stokes I':
            im = ax.contourf(ra_principle, dec_principle, self.I, cmap=color, levels=50)
        elif value_name == 'Stokes Q': # maybe mistmatch
            im = ax.contourf(ra_principle, dec_principle, q_principle, cmap=color, levels=50)
            contours1 = ax.contour(ra_principle, dec_principle, q_principle, colors='red', linewidths=.5, levels=[3 * self.noise_Q])
            contours2 = ax.contour(ra_principle, dec_principle, q_principle, colors='blue', linewidths=.5, levels=[-3 * self.noise_Q])
        elif value_name == 'Stokes U':
            im = ax.contourf(ra_principle, dec_principle, u_principle, cmap=color, levels=50)
            contours1 = ax.contour(ra_principle, dec_principle, u_principle, colors='red', linewidths=.5, levels=[3 * self.noise_U])
            contours2 = ax.contour(ra_principle, dec_principle, u_principle, colors='blue', linewidths=.5, levels=[-3 * self.noise_U])

        ax.set_xlabel(r'$\Delta$RA (arcsec)')
        ax.set_ylabel(r'$\Delta$Dec (arcsec)')
        ax.set_aspect('equal')

        #plot axis
        slope = np.tan(np.radians(90 - self.angle_pa))
        ax.axline((0, 0), slope=0, linestyle='--', linewidth=0.5, color='black', alpha=.5)
        ax.axline((0, 0), slope=np.inf, linestyle='--',linewidth=0.5, color='black',alpha=.5)


        ax.contour(ra_principle,dec_principle , self.I, colors='white', linewidths=.4, alpha=.4,
                   levels=[3 * self.noise_I, 10 * self.noise_I, 25 * self.noise_I, 50 * self.noise_I,
                           100 * self.noise_I,
                           200 * self.noise_I,
                           325 * self.noise_I, 500 * self.noise_I, 1000 * self.noise_I])


        #uncomment to ensure principle frame transformation is correct
        #ax.plot(x, maj, color='black', linestyle='-', linewidth=2)


        #plot au bar
        ax.quiver(1,-self.window+1, 0.82644628099,0,  pivot='middle', angles='xy', color='black', scale=1, scale_units='xy', headwidth=1e-10,
                  headlength=1e-10, headaxislength=1e-10, width=0.005)
        ax.annotate(text='100 AU', xy=(1,1.8), xytext=(1.4,-self.window+1.1), textcoords='data', fontsize=7, color='black', weight='bold')

        #plot beam
        self.plot_beam(ax,self.bmaj, self.bmin, self.bpa - (self.angle_pa -self.bpa))

        ax.set_title(f'{value_name} {self.band} In Principle Frame')

        cbar = plt.colorbar(im)
        cbar.set_label(f'{value_name} [{unit}]')

        ax.set_xlim(self.window, -self.window)
        ax.set_ylim(-self.window, self.window)
        return ax

    def plot_vect_radius(self,ax, comparison=False, principle_frame = False):
        #sample points along radius
        ra_sample, dec_sample, vra_azi, vdec_azi = self.create_radius(principle_frame)

        if principle_frame:
            eta = np.radians(self.angle_pa - 90)  # consistent with principle_frame()


            #   (inverse of the principle_frame() transformation)
            ra_sky = ra_sample * np.cos(eta) + dec_sample * np.sin(eta)
            dec_sky = -ra_sample * np.sin(eta) + dec_sample * np.cos(eta)
            points = np.column_stack((dec_sky, ra_sky))


            Q_interp = rgi((self.dec, self.ra), self.Q_princ,
                           bounds_error=False, fill_value=np.nan)
            U_interp = rgi((self.dec, self.ra), self.U_princ,
                           bounds_error=False, fill_value=np.nan)
            I_interp = rgi((self.dec, self.ra), self.I,
                           bounds_error=False, fill_value=np.nan)

            Q_sample = Q_interp(points)
            U_sample = U_interp(points)
            I_sample = I_interp(points)

            p_sample, pf_sample, chi_error = self.calculate_P_debiased(
                Q_sample, U_sample, I_sample)



        else:
            #instantiate rectangular grid interpolator for normal coords
            Q_interp = rgi((self.dec,self.ra), self.Q, bounds_error=False, fill_value=np.nan) #(dec, ra)
            U_interp = rgi((self.dec,self.ra), self.U, bounds_error=False, fill_value=np.nan)
            I_interp = rgi((self.dec,self.ra), self.I, bounds_error=False, fill_value=np.nan)

            points = np.column_stack((dec_sample, ra_sample))
                    # interpolate from given points
            Q_sample = Q_interp(points)
            U_sample = U_interp(points)
            I_sample = I_interp(points)

            Q_sample = Q_sample.reshape(dec_sample.shape)
            U_sample = U_sample.reshape(dec_sample.shape)
            I_sample = I_sample.reshape(dec_sample.shape)

            p_sample, pf_sample, chi_error = self.calculate_P_debiased(Q_sample, U_sample, I_sample)
            print(f'sky frame p: {p_sample, pf_sample}')
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
        return pa_arr, comparison_angle, chi_error

    def calculate_P_debiased(self,q_sample, u_sample, I_sample):
        P_m = np.sqrt(q_sample ** 2 + u_sample ** 2)
        P_debiased = np.full_like(q_sample, np.nan)

        P_simple = np.sqrt(np.maximum(q_sample ** 2 + u_sample ** 2 - self.noise_P ** 2, 0))

        pm_values = np.linspace(0, 5 * self.noise_P, 50)
        p_values = np.linspace(0, self.I.max(), 100)
        results = np.full_like(pm_values, np.nan)

        for i in range(len(pm_values)):
            results[i] = np.argmax(
                p_values / self.noise_P ** 2
                * i0((p_values * pm_values[i]) / (self.noise_P ** 2))
                * np.exp(-(pm_values[i] ** 2 + p_values ** 2) / (2 * self.noise_P ** 2))
            )

        interp_obj = make_interp_spline(pm_values, results)
        P_debiased = np.where(P_m < 5 * self.noise_P, interp_obj(P_m), P_simple)

        mask = (
                (I_sample > 3 * self.noise_I) &
                (P_debiased > 3 * self.noise_P)
        )
        pf_debiased = np.full_like(I_sample, np.nan)
        # pf = P_simple/I_arr
        pf_debiased = np.where(mask, P_debiased/I_sample, np.nan)
        chi_error = np.where(mask, 28.65 * self.noise_P/P_debiased, np.nan)

        return P_debiased, pf_debiased, chi_error

    def create_radius(self , principle=False):

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
        ra_im = ra_disk
        dec_im = dec_disk * np.cos(np.radians(self.angle_incl))

        vra_im =  ra_disk_azi
        vdec_im = dec_disk_azi * np.cos(np.radians(self.angle_incl))

        if not principle:
            #rotate
            ra_im0 = ra_im
            dec_im0 = dec_im
            vra_im0 = vra_im
            vdec_im0 = vdec_im

            eta = -np.radians(self.angle_pa-90)
            #PROBLEM IS PROBABLY HERE WHERE AXIS SWITCH IS MISMATCHING COORDS FROM PRINCIPLE FRAME TO SKY
            # x = x0 * np.cos(angle_pa) - np.sin(angle_pa) * y0
            # y = y0 * np.cos(angle_pa) + np.sin(angle_pa) * x0

            # ra_im = ra_im0 * np.cos(eta) - np.sin(eta) * dec_im0
            # dec_im = dec_im0 * np.cos(eta) + np.sin(eta) * ra_im
            #
            # vdec_im = vdec_im0 * np.cos(eta) + np.sin(eta) * ra_im0
            #vra_im = vra_im0 * np.cos(eta) - np.sin(eta) * vdec_im
            dec_im = ra_im0 * np.sin(eta) + np.cos(eta) * dec_im0
            ra_im = -dec_im0 * np.sin(eta) + np.cos(eta) * ra_im0

            vdec_im = vra_im0 * np.sin(eta) + np.cos(eta) * vdec_im0
            vra_im = -vdec_im0 * np.sin(eta) + np.cos(eta) * vra_im0
            #plt.plot(ra_im, dec_im, color='red')
        return ra_im, dec_im, vra_im, vdec_im


    def compare_angles(self, obs_angle, exp_angle, chi_error):
        diff_angle = np.degrees(obs_angle - exp_angle)
        diff_angle = (diff_angle + 90) % 180 - 90
        ax = plt.gca()
        ax.hist(diff_angle, alpha = .7, color = 'blue')
        ax.set_xlabel('Difference of Observed and Expected Angle')
        ax.set_ylabel('Frequency')

        mask = np.isfinite(diff_angle)
        data_mean = np.mean(diff_angle[mask])
        typical_error = np.median(chi_error[mask])
        print(np.sqrt(np.mean(data_mean ** 2)))

        ax.axvline(float(data_mean), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {data_mean:.2f}')
        ax.axvline(float(3*typical_error), color='gray', linestyle='dashed', linewidth=2, )
        ax.axvline(float(0 - 3*typical_error), color='gray', linestyle='dashed', linewidth=2,)

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


    def principle_frame(self):
        #debugging only
        x = np.linspace(-4, 4, 100)
        slope = np.tan(np.radians(90 - self.angle_pa))
        maj = slope * x
        ##


        angle_pa = np.radians(self.angle_pa-90)

        #x is plotted vertically
        x_principle = self.ra * np.sin(angle_pa) + np.cos(angle_pa) * self.dec
        y_principle = -self.dec * np.sin(angle_pa) + np.cos(angle_pa) * self.ra

        maj_principle = maj * np.sin(angle_pa) - np.cos(angle_pa) * x
        x_principle_test = x * np.sin(angle_pa) + np.cos(angle_pa) * maj

        #take note the grid is never rotated with, so q' and u' are left on the original grid
        q_principle = self.U * np.sin(2*angle_pa) + np.cos(2*angle_pa) * self.Q
        u_principle = -self.Q * np.sin(2*angle_pa) + np.cos(2*angle_pa) * self.U

        self.set_principle_frame( q_principle, u_principle, x_principle, y_principle)

        return x_principle, y_principle, q_principle, u_principle, maj_principle, x_principle_test

    def test1(self, obs_princ, obs_sky):
        eta = -np.radians(90-self.angle_pa)
        pa_princ_back = obs_princ + eta
        diff = np.degrees(obs_sky - pa_princ_back)
        diff = (diff + 90) % 180 - 90
        print(np.nanmean(diff))
        print(np.nanstd(diff))
        plt.hist(diff[np.isfinite(diff)], bins=50)
        plt.show()

    def test2(self):
        eta = -np.radians(90 - self.angle_pa)
        pa_sky = 0.5 * np.arctan2(self.U, self.Q)

        pa_from_princ = (
                0.5 * np.arctan2(self.U_princ, self.Q_princ)
                + eta
        )

        diff = np.degrees(pa_sky - pa_from_princ)
        diff = (diff + 90) % 180 - 90

        print(np.nanmean(diff))
        print(np.nanstd(diff))

    def test3(self):
        self.principle_frame()
        ra_sky, dec_sky, _, _ = self.create_radius(False)


        Q_interp = rgi((self.dec, self.ra), self.Q,
                       bounds_error=False, fill_value=np.nan)

        U_interp = rgi((self.dec, self.ra), self.U,
                       bounds_error=False, fill_value=np.nan)

        points = np.column_stack((dec_sky, ra_sky))

        Q_sky = Q_interp(points)
        U_sky = U_interp(points)

        pa_sky = 0.5 * np.arctan2(U_sky, Q_sky)

        x_princ, y_princ, _, _ = self.create_radius(True)

        Q_interp = rgi((self.dec, self.ra), self.Q_princ,
                       bounds_error=False, fill_value=np.nan)

        U_interp = rgi((self.dec, self.ra), self.U_princ,
                       bounds_error=False, fill_value=np.nan)
        eta = np.radians(self.angle_pa - 90)
        ra_sky = x_princ * np.cos(eta) + y_princ * np.sin(eta)
        dec_sky = -x_princ * np.sin(eta) + y_princ * np.cos(eta)
        points = np.column_stack((dec_sky, ra_sky))



        Q_princ = Q_interp(points)
        U_princ = U_interp(points)

        pa_princ = 0.5 * np.arctan2(U_princ, Q_princ)


        pa_princ_back = pa_princ + eta
        diff_obs = np.degrees(pa_sky - pa_princ_back)
        diff_obs = (diff_obs + 90) % 180 - 90

        print("mean =", np.nanmean(diff_obs))
        print("std  =", np.nanstd(diff_obs))

        plt.hist(diff_obs[np.isfinite(diff_obs)], bins=50)
        plt.title("Observed angle field invariance")
        plt.show()

