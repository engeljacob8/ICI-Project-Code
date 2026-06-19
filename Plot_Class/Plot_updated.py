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
from scipy.optimize import minimize
from scipy import stats

class Plot_updated:
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
        self.set_principle_frame(None,None,None, None, None)
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
        self.noise_P = np.sqrt(max((noise_Q**2 + noise_U**2)/2,0))

    def set_stokes(self, I, Q, U):
        self.I = np.array(I)
        self.Q = np.array(Q)
        self.U = np.array(U)

    def set_principle_frame(self,  Q_princ, U_princ, I_princ, x_princ, y_princ ):

        self.Q_princ = Q_princ
        self.U_princ = U_princ
        self.I_princ = I_princ
        self.x_princ = x_princ
        self.y_princ = y_princ




    def plot_image (self, value, value_name, unit,overlay_extent=None,overlay=False,
              contour=False, ax=None,color='viridis'):
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




        slope = np.tan(np.radians(90-self.angle_pa))
        ax.axline((0,0),slope=slope, linestyle='--' ,linewidth=0.5, color='black',alpha=.5)
        ax.axline((0, 0), slope= -1/slope, linestyle='--',linewidth=0.5, color='black',alpha=.5)

        if overlay:
            ax.set_title(f'Polarization at band {self.band} + 1.2mm Dust')
        else:
            ax.set_title(f'{value_name} {self.band}')

        self.plot_beam(ax)

        ax.quiver(1,-self.window+1, 0.82644628099,0,  pivot='middle', angles='xy', color='black', scale=1, scale_units='xy', headwidth=1e-10,
                  headlength=1e-10, headaxislength=1e-10, width=0.005)
        ax.annotate(text='100 AU', xy=(1,1.8), xytext=(1.4,-self.window+1.1), textcoords='data', fontsize=7, color='black', weight='bold')


        contours = ax.contour(self.ra, self.dec, self.I, colors='white', linewidths=.6, alpha=.5,
                                  levels=[3 * self.noise_I, 10 * self.noise_I, 25 * self.noise_I, 50 * self.noise_I, 100 * self.noise_I,
                                          200 * self.noise_I,
                                          325 * self.noise_I, 500 * self.noise_I, 1000 * self.noise_I])


        cbar = plt.colorbar(im)
        cbar.set_label(f'{value_name} [{unit}]')

        ax.set_xlim(self.window, -self.window)
        ax.set_ylim(-self.window, self.window)
        return ax

    def plot_principle_frame(self, value_name, unit, ax=None, color='viridis' ):
        if ax is None:
            ax = plt.gca()
        q_principle, u_principle, I_principle= self.principle_frame()

        if value_name == 'Stokes I':
            im = ax.contourf(self.ra, self.dec, I_principle, cmap=color, levels=50)
        elif value_name == 'Stokes Q':
            im = ax.contourf(self.ra, self.dec, q_principle, cmap=color, levels=50)
            contours1 = ax.contour(self.ra, self.dec, q_principle, colors='red', linewidths=.5, levels=[3 * self.noise_Q])
            contours2 = ax.contour(self.ra, self.dec, q_principle, colors='blue', linewidths=.5, levels=[-3 * self.noise_Q])
        elif value_name == 'Stokes U':
            im = ax.contourf(self.ra, self.dec, u_principle, cmap=color, levels=50)
            contours1 = ax.contour(self.ra, self.dec, u_principle, colors='red', linewidths=.5, levels=[3 * self.noise_U])
            contours2 = ax.contour(self.ra, self.dec, u_principle, colors='blue', linewidths=.5, levels=[-3 * self.noise_U])

        ax.set_xlabel(r'$\Delta$RA (arcsec)')
        ax.set_ylabel(r'$\Delta$Dec (arcsec)')
        ax.set_aspect('equal')

        #plot axis

        ax.axline((0, 0), slope=0, linestyle='--', linewidth=0.5, color='black', alpha=.5)
        ax.axline((0, 0), slope=np.inf, linestyle='--',linewidth=0.5, color='black',alpha=.5)


        ax.contour(self.ra,self.dec , I_principle, colors='white', linewidths=.4, alpha=.4,
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


        ax.set_title(f'{value_name} {self.band} In Principal Frame')

        cbar = plt.colorbar(im)
        cbar.set_label(f'{value_name} [{unit}]')

        ax.set_xlim(self.window, -self.window)
        ax.set_ylim(-self.window, self.window)
        return ax

    def plot_vectors(self, ax, comparison=False, principle_frame=False, sample= True):
        # sample points along radius
        ra_sample, dec_sample, vra_azi, vdec_azi = self.create_radius(principle_frame)

        if principle_frame:
            eta = np.radians(90-self.angle_pa )
            # dec_im = ra_im0 * np.sin(eta) + np.cos(eta) * dec_im0
            # ra_im = -dec_im0 * np.sin(eta) + np.cos(eta) * ra_im0

            dec_sky = ra_sample * np.sin(eta) + dec_sample * np.cos(eta)
            ra_sky = ra_sample * np.cos(eta) - dec_sample * np.sin(eta)

            points = np.column_stack((dec_sky, ra_sky))
            # points = np.column_stack((dec_sample, ra_sample))

            Q_interp = rgi((self.dec, self.ra), self.Q_princ,
                           bounds_error=False, fill_value=np.nan)
            U_interp = rgi((self.dec, self.ra), self.U_princ,
                           bounds_error=False, fill_value=np.nan)
            I_interp = rgi((self.dec, self.ra), self.I,
                           bounds_error=False, fill_value=np.nan)

            Q_sample = Q_interp(points).reshape(dec_sky.shape)
            U_sample = U_interp(points).reshape(dec_sky.shape)
            I_sample = I_interp(points).reshape(dec_sky.shape)

            p_sample, pf_sample, chi_error = self.calculate_P_debiased(Q_sample, U_sample, I_sample)
        else:
            # instantiate rectangular grid interpolator for normal coords
            Q_interp = rgi((self.dec, self.ra), self.Q, bounds_error=False, fill_value=np.nan)  # (dec, ra)
            U_interp = rgi((self.dec, self.ra), self.U, bounds_error=False, fill_value=np.nan)
            I_interp = rgi((self.dec, self.ra), self.I, bounds_error=False, fill_value=np.nan)

            points = np.column_stack((dec_sample, ra_sample))
            # interpolate from given points
            Q_sample = Q_interp(points)
            U_sample = U_interp(points)
            I_sample = I_interp(points)

            Q_sample = Q_sample.reshape(dec_sample.shape)
            U_sample = U_sample.reshape(dec_sample.shape)
            I_sample = I_sample.reshape(dec_sample.shape)

            p_sample, pf_sample, chi_error = self.calculate_P_debiased(Q_sample, U_sample, I_sample)

        # calculate polarization angle

        pa_arr = .5 * np.arctan2(U_sample, Q_sample)  # this angle is east of north

        # cartesian
        # arcsec for 1% polarization
        scale = .5  # arcsec/percent
        d_ra = scale * pf_sample * 100 * np.sin(pa_arr)  # dec on x ax
        d_dec = scale * pf_sample * 100 * np.cos(pa_arr)
        pa_arr = np.arctan2(d_ra, d_dec)

        comparison_angle = np.arctan2(vra_azi, vdec_azi)
        vra = np.sin(comparison_angle)
        vdec = np.cos(comparison_angle)




        q = ax.quiver(ra_sample, dec_sample, d_ra, d_dec, pivot='middle', angles='xy', color='red', scale=1,
                      scale_units='xy', headwidth=1e-10,
                      headlength=1e-10, headaxislength=1e-10, width=0.005)
        ax.quiverkey(q, X=self.window-1, Y=self.window-.5, U=scale, label=f'{scale} % Polarization', labelpos='E',coordinates= 'data')

        if comparison:
            ax.quiver(ra_sample, dec_sample, vra, vdec, pivot='middle', angles='xy', color='white', scale=5,
                      scale_units='xy', headwidth=1e-10,
                      headlength=1e-10, headaxislength=1e-10, width=0.005)
            ax.set_title('Expected vs Observed Polarization')

        return pa_arr, comparison_angle, chi_error


    def create_radius(self , principle=False, sample=False):
        if sample:
            radii = np.array([.75])
        else:
            radii = np.array([.25* self.bmaj, .75 * self.bmaj, 1.25 * self.bmaj, 1.75 * self.bmaj])
        arc_spacing = self.bmaj / 2

        ra_list = []
        dec_list = []

        for r in radii:
            #approximation for circular distances of 1/2 bmaj - can be improved with elliptical circum
            spacing = max(4, int(np.ceil(2 * np.pi * r / arc_spacing)))
            phi = np.linspace(-np.pi/2, 3* np.pi/2, spacing, endpoint=False)

            ra_list.append(r * np.sin(phi))
            dec_list.append(r * np.cos(phi))

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
        if sample:
            return ra_im, dec_im, phi

        if not principle:
            #rotate
            ra_im0 = ra_im
            dec_im0 = dec_im
            vra_im0 = vra_im
            vdec_im0 = vdec_im

            #angle is north of east
            eta = np.radians(90-self.angle_pa)

            dec_im = ra_im0 * np.sin(eta) + np.cos(eta) * dec_im0
            ra_im = -dec_im0 * np.sin(eta) + np.cos(eta) * ra_im0

            vdec_im = vra_im0 * np.sin(eta) + np.cos(eta) * vdec_im0
            vra_im = -vdec_im0 * np.sin(eta) + np.cos(eta) * vra_im0
            #plt.plot(ra_im, dec_im, color='red')
        return ra_im, dec_im, vra_im, vdec_im



    def principle_frame(self):
        eta = -np.radians(90-self.angle_pa)

        #take note the grid is never rotated with, so q' and u' are left on the original grid
        q_principle = self.U * np.sin(2*eta) + np.cos(2*eta) * self.Q
        u_principle = -self.Q * np.sin(2*eta) + np.cos(2*eta) * self.U

        x_principle = self.ra * np.sin(eta) + np.cos(eta) * self.dec
        y_principle = -self.dec * np.sin(eta) + np.cos(eta) * self.ra

        # pf and I principle DO live on the principle frame grid
        #change angle_pa
        eta = -eta
        interp_I = rgi((self.dec, self.ra), self.I, bounds_error=False, fill_value=np.nan )
        RAp, DECp = np.meshgrid( self.ra, self.dec )


        #possible point of error
        DEC = RAp * np.sin(eta) + np.cos(eta) * DECp
        RA = -DECp * np.sin(eta) + np.cos(eta) * RAp

        points = np.column_stack([ DEC.ravel(), RA.ravel() ])
        I_principle = interp_I(points).reshape(RAp.shape)

        self.set_principle_frame(q_principle, u_principle,I_principle, x_principle, y_principle)

        return  q_principle, u_principle,I_principle



    def sample_azimuth(self):
        self.principle_frame()
        ra_sample, dec_sample, phi = self.create_radius(principle=True, sample=True)

        eta = np.radians(90 - self.angle_pa)

        dec_sky = ra_sample * np.sin(eta) + dec_sample * np.cos(eta)
        ra_sky = ra_sample * np.cos(eta) - dec_sample * np.sin(eta)

        points = np.column_stack((dec_sky, ra_sky))

        Q_interp = rgi((self.dec, self.ra), self.Q_princ,
                       bounds_error=False, fill_value=np.nan)
        U_interp = rgi((self.dec, self.ra), self.U_princ,
                       bounds_error=False, fill_value=np.nan)
        I_interp = rgi((self.dec, self.ra), self.I,
                       bounds_error=False, fill_value=np.nan)

        Q_sample = Q_interp(points).reshape(dec_sky.shape)
        U_sample = U_interp(points).reshape(dec_sky.shape)
        I_sample = I_interp(points).reshape(dec_sky.shape)

        p_sample, pf_sample, chi_error = self.calculate_P_debiased(Q_sample, U_sample, I_sample)
        ax = plt.gca()

        sigma_q = np.sqrt(
            (self.noise_Q / I_sample) ** 2 +
            (Q_sample * self.noise_I / I_sample ** 2) ** 2
        )
        sigma_u = np.sqrt(
            (self.noise_U / I_sample) ** 2 +
            (U_sample * self.noise_I / I_sample ** 2) ** 2
        )
         #measured east of north
        q_prime = np.where(pf_sample > 0, Q_sample / I_sample,np.nan)
        u_prime = np.where(pf_sample > 0, U_sample / I_sample, np.nan)
        mask_q = np.isfinite(q_prime) & np.isfinite(sigma_q)

        phi_q = phi[mask_q]
        q_prime = q_prime[mask_q]
        sigma_q = sigma_q[mask_q]

        mask_u = np.isfinite(u_prime) & np.isfinite(sigma_u)

        phi_u = phi[mask_u]
        u_prime = u_prime[mask_u]
        sigma_u = sigma_u[mask_u]


        ax.errorbar(
            np.degrees(phi_q),
            q_prime *100,
            yerr=sigma_q*100,
            fmt='o',
            color='black',
            alpha=0.7
        )
        self.plot_model_q(ax, phi_q, q_prime, sigma_q)

        ax.set_xlabel(r' $\phi$ (deg)', fontsize=12)
        ax.set_ylabel(r"q' [%]", fontsize=12)

        ax.set_xlim(-95, 270)
        ax.set_xticks(np.arange(-90, 271, 45))


        ax.grid(True, alpha=0.3)

        ax.set_title(r"")

        plt.tight_layout()
        plt.show()

        ax1 = plt.gca()
        ax1.errorbar(
            np.degrees(phi_u),
            u_prime * 100,
            yerr=sigma_u*100,
            fmt='o',
            color='black',
            alpha=0.7
        )
        self.plot_model_u(ax1, phi_u, u_prime, sigma_u)

        ax1.set_xlabel(r' $\phi$ (deg)', fontsize=12)
        ax1.set_ylabel(r"u' [%]", fontsize=12)

        ax1.set_xlim(-95, 270)

        ax1.set_xticks(np.arange(-90, 271, 45))

        ax1.grid(True, alpha=0.3)

        ax1.set_title(r"")

        plt.tight_layout()
        plt.show()

    def plot_model_q(self,ax, phi, q_sample, sigma_q):
        theta =  np.linspace(-np.pi/2, 3* np.pi/2, 100)
        ax.axline([0,0], slope=0, color='black', linestyle='--')

        t, s = self.mle_estimation_q(phi, q_sample, sigma_q)
        q_model = s + t * (np.cos(np.radians(self.angle_incl))**2 * np.sin(theta)**2 - np.cos(theta)**2)
        print(t* 100)
        print(s* 100)
        ax.plot(np.degrees(theta),q_model * 100, '-', color='green')

        ax.axline([0,s*100], slope=0, color='green', linestyle='--')



    def plot_model_u(self, ax, phi, u_sample, sigma_u):
        theta = np.linspace(-np.pi/2, 3* np.pi/2, 100)
        ax.axline([0,0], slope=0, color='black', linestyle='--')

        t = self.mle_estimation_u(phi, u_sample, sigma_u)
        print(t * 100)

        u_model = -t * np.cos(np.radians(self.angle_incl)) * np.sin(2* theta)
        ax.plot(np.degrees(theta), u_model * 100, '-', color='green')


    def mle_estimation_q(self, phi, q_sample, sigma_q):
        result = minimize( lambda x: -self.likelihood_function_q(x,phi,q_sample,sigma_q), x0= [.1, 0.0])
        t , s = result.x
        q_model = s + t * (np.cos(np.radians(self.angle_incl)) ** 2 * np.sin(phi) ** 2 - np.cos(phi) ** 2)

        chi2 = np.sum(((q_sample - q_model) / sigma_q) ** 2)
        p_value = stats.chi2.sf(chi2, len(q_sample))
        print(f' p_value of reduced chi^2 test for q model {p_value} ')
        return t, s

    def mle_estimation_u(self, phi, u_sample, sigma_u):
        result = minimize(lambda x: -self.likelihood_function_u(x, phi, u_sample, sigma_u), x0=[.1])
        t = result.x[0]

        u_model = -t * np.cos(np.radians(self.angle_incl)) * np.sin(2 * phi)

        chi2 = np.sum(((u_sample-u_model)/sigma_u)**2)
        p_value = stats.chi2.sf(chi2, len(u_sample))
        print(f' p_value of reduced chi^2 test for u model {p_value} ')
        return t


    def likelihood_function_u(self,params, phi, u_sample, sigma_u):
        t = params[0]
        u_model = -t * np.cos(np.radians(self.angle_incl)) * np.sin(2 * phi)

        return -.5 * np.sum(((u_sample-u_model)**2 / sigma_u**2) + np.log(2 * np.pi * sigma_u**2))

    def likelihood_function_q(self,params, phi, q_sample, sigma_q):
        t, s = params
        q_model = s + t * (np.cos(np.radians(self.angle_incl)) ** 2 * np.sin(phi) ** 2 - np.cos(phi) ** 2)

        return -.5 * np.sum(((q_sample-q_model)**2 / sigma_q**2) + np.log(2 * np.pi * sigma_q**2))
















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

    def compare_angles(self, obs_angle, exp_angle, chi_error):

        diff_angle = np.degrees(obs_angle - exp_angle)
        diff_angle = (diff_angle + 90) % 180 - 90
        ax = plt.gca()
        ax.hist(diff_angle, alpha = .7, color = 'blue')
        ax.set_xlabel('Difference of Observed and Expected Angle')
        ax.set_ylabel('Frequency')

        diff2 = np.degrees(
            0.5 * np.arctan2(
                np.sin(2 * (obs_angle - exp_angle)),
                np.cos(2 * (obs_angle - exp_angle))
            )
        )
        mask = np.isfinite(diff_angle)
        data_mean = np.mean(diff_angle[mask])
        typical_error = np.median(chi_error[mask])
        print(np.sqrt(np.mean(diff_angle[mask] ** 2)))

        ax.axvline(float(data_mean), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {data_mean:.2f}')
        ax.axvline(float(3*typical_error), color='gray', linestyle='dashed', linewidth=2, )
        ax.axvline(float(0 - 3*typical_error), color='gray', linestyle='dashed', linewidth=2,)

        plt.show()
        ax2 = plt.gca()
        ax2.hist(diff2, alpha=.7, color='blue')
        ax2.set_xlabel('Difference of Observed and Expected Angle')
        ax2.set_ylabel('Frequency')
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

    def test1(self, obs_princ, obs_sky):
        eta = np.radians(90-self.angle_pa)
        pa_princ_back = obs_princ + eta
        diff = np.degrees(obs_sky - pa_princ_back)
        diff = (diff + 90) % 180 - 90
        print(np.nanmean(diff))
        print(np.nanstd(diff))


    def test2(self):
        eta = np.radians(90 - self.angle_pa)
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

        plt.hist(diff_obs[np.isfinite(diff_obs)], bins=10)
        plt.title("Observed angle field invariance")
        plt.show()