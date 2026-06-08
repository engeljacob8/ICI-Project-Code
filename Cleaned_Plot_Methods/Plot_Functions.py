import numpy as np
import matplotlib.pyplot as plt
from radio_beam import Beam
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator as rgi


def image_plot_pf(pf,extent, **kwargs ):
    ax = plt.gca()
    im = ax.imshow( pf, cmap='viridis', extent= extent  )
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Polarization Fraction [%]')

    if 'bmaj' in kwargs:
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
         ax.add_artist(artist)

    #plot_radius(34.9, ax)

    plt.xlim(3, -3) #manually plots x in other direction
    plt.ylim(-3, 3)
    plt.show()

def cont_plot_pf(pf,  ra, dec, **kwargs):
    ax = plt.gca()
    im = ax.contourf(ra, dec ,pf, cmap='viridis',  levels = 50 ) #RA plotted on x-axis
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Polarization Fraction [%]')

    if 'bmaj' in kwargs:
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
         ax.add_artist(artist)

    plt.xlim(3, -3) #manually plots x in other direction
    plt.ylim(-3, 3)
    plt.show()

def im_plot_p(P, extent, noise_I, I, **kwargs):
    ax = plt.gca()
    im = ax.imshow( P,cmap='viridis' , extent= extent)
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours = ax.contour(I, colors='white', linewidths=.5, extent=extent ,     #plots I contour levels
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000 * noise_I])
    ax.clabel(contours, inline=1, fontsize=10)

    if 'bmaj' in kwargs:
        artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
        ax.add_artist(artist)

    cbar = plt.colorbar(im)
    cbar.set_label('Linear Polarization Intensity [Jy/beam]')

    plt.xlim(3,-3)
    plt.ylim(-3,3)
    plt.show()

def cont_plot_p(P, I, noise_I, ra, dec, **kwargs):
    ax = plt.gca()
    im = ax.contourf(ra, dec, P, cmap='viridis', levels=50)
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours = ax.contour(ra, dec, I, colors='white', linewidths=.5,
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000 * noise_I])
    ax.clabel(contours, inline=1, fontsize=10)

    if 'bmaj' in kwargs:
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
         ax.add_artist(artist)

    plot_radius(34.9, 85.8, ax)

    cbar = plt.colorbar(im)
    cbar.set_label('P [Jy/beam]')

    plt.xlim(3,-3)
    plt.ylim(-3, 3)
    plt.show()

def plot_q(Q, ra, dec, I, noise_I,**kwargs):
    ax = plt.gca()
    im = ax.contourf(ra,dec, Q,cmap='viridis', levels=50, )
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours = ax.contour(ra, dec, I, colors='white', linewidths=.5,
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000 * noise_I])
    ax.clabel(contours, inline=1, fontsize=10)

    if 'bmaj' in kwargs:
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
         ax.add_artist(artist)

    cbar = plt.colorbar(im)
    cbar.set_label('Stokes Q [Jy/beam]')

    plt.xlim(3,-3) # manually swap
    plt.ylim(-3,3)
    plt.show()

def plot_i(I, ra, dec, noise_I,**kwargs):
    ax = plt.gca()
    im = ax.contourf(ra, dec,I,cmap='viridis', levels=50, )
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours = ax.contour(ra, dec, I, colors= 'white', linewidths= 1.5,
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000*noise_I])
    ax.clabel(contours, inline=1, fontsize=10)

    cbar = plt.colorbar(im)

    cbar.set_label('Stokes I [Jy/beam]')
    if 'bmaj' in kwargs:
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'], ax)


    #if 'U_arr' in kwargs:
        #plot_vectors(ra, dec, ax, kwargs['bmaj'], kwargs['bmin'], kwargs['U_arr'], kwargs['Q_arr'], kwargs['pf'])

    #plot_radius(34.9, 85.8, ax)
    plot_vect_on_radius(ra, dec, ax, kwargs['bmaj'], kwargs['bmin'], kwargs['U_arr'], kwargs['Q_arr'], kwargs['pf'])

    plt.xlim(3,-3) #manually swapped direction
    plt.ylim(-3,3)
    plt.show()


# def plot_i_test(I, ra, dec, noise_I,**kwargs):
#     ax = plt.gca()
#     im = ax.contourf(ra, dec,I,cmap='viridis', levels=50, )
#     ax.set_xlabel(r'$\Delta$RA (arcsec)')
#     ax.set_ylabel(r'$\Delta$Dec (arcsec)')
#     ax.set_aspect('equal')
#
#     contours = ax.contour(ra, dec, I, colors= 'white', linewidths= 1.5,
#                           levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
#                                   325 * noise_I, 500 * noise_I, 1000*noise_I])
#     ax.clabel(contours, inline=1, fontsize=10)
#
#     cbar = plt.colorbar(im)
#
#     cbar.set_label('Stokes I [Jy/beam]')
#     if 'bmaj' in kwargs:
#          artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
#          ax.add_artist(artist)
#
#     if 'Ux' in kwargs:
#         ax.quiver(ra, dec, kwargs['Ux'], kwargs['Uy'],pivot='middle', scale=10, units='xy' )
#
#     plt.xlim(3,-3) #manually swapped direction
#     plt.ylim(-3,3)
#     plt.show()

def plot_u(U, ra, dec, I, noise_I, **kwargs):
    ax = plt.gca()
    im = ax.contourf(ra, dec, U, cmap='viridis', levels=50, )
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours1 = ax.contour(ra, dec, U, colors= 'red', linewidths=.5, levels = [3*kwargs['noise_U']])
    contours2 = ax.contour(ra, dec, U, colors='blue', linewidths=.5, levels=[3 * -kwargs['noise_U']])

    contours = ax.contour(ra, dec, I, colors='white', linewidths=1.5,
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000 * noise_I])

    ax.clabel(contours,inline= 1,  fontsize=10, colors='red') # not working lol


    if 'bmaj' in kwargs:
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
         ax.add_artist(artist)

    cbar = plt.colorbar(im)
    cbar.set_label('Stokes U [Jy/beam]')

    plt.xlim(3, -3) # manually swapped direction
    plt.ylim(-3, 3)
    plt.show()

def plot_beam(bmaj, bmin,bpa,ax):
    my_beam = Beam( bmaj * u.arcsec, bmin*u.arcsec, bpa*u.degree)
    y_cen_pix, x_cen_pix = -2,2
    pixscale = u.arcsec
    ellipse_artist = my_beam.ellipse_to_plot(x_cen_pix, y_cen_pix, pixscale, facecolor='black')
    ax.add_artist(ellipse_artist)


def plot_vectors(x_axis, y_axis,ax, bmaj, bmin , U_arr, Q_arr, pf):
    spacing_y = bmin/2
    spacing_x = bmaj/2
    #sample points at 1/2 beams
    y_sample = np.arange(y_axis.min(), y_axis.max(), spacing_y)
    x_sample = np.arange(x_axis.min(), x_axis.max(), spacing_x)
    y_points, x_points = np.meshgrid(y_sample, x_sample) #done opposite to match the fits data



    #instantiate rectangular grid interpolator
    Q_interp = rgi((y_axis, x_axis), Q_arr, bounds_error=False, fill_value=np.nan)
    U_interp = rgi( (y_axis, x_axis), U_arr, bounds_error=False, fill_value=np.nan)
    pf_interp = rgi( (y_axis, x_axis), pf, bounds_error=False, fill_value=np.nan)

    points = np.column_stack((y_points.flatten(), x_points.flatten()))
    #interpolate from given points
    Q_sample = Q_interp(points)
    U_sample = U_interp(points)
    pf_sample = pf_interp(points)

    Q_sample = Q_sample.reshape(y_points.shape)
    U_sample = U_sample.reshape(y_points.shape)
    pf_sample = pf_sample.reshape(y_points.shape)

    #calculate polarization angle
    pa_arr = .5 * np.arctan2(U_sample, Q_sample)  # this angle is east of north
    # cartesian
    # .3 arcsec for 1% polarization
    d_x = .3 * pf_sample*100 * -np.sin(pa_arr)
    d_y = .3 * pf_sample*100 * np.cos(pa_arr)

    #yy plotted opposite since ra is plotted on x-axis in conventional graph
    ax.quiver( x_points,y_points, d_x,d_y, pivot='middle',color='red', scale=1, scale_units='xy',headwidth=1e-10, headlength=1e-10, headaxislength=1e-10, width=0.005)


def plot_radius(angle_incl, angle_pa ,ax):
    phi = np.linspace(0, 2*np.pi, 100)
    plots = 4
    for i in range(plots):
        radius_plots = i # arc sec for now

        x0 = radius_plots * np.cos(phi)
        y0= radius_plots * np.cos(np.radians(angle_incl)) * np.sin(phi)

        # rotate - works for when east is to the left
        angle_pa = np.radians(angle_pa)
        x = x0 * np.cos(angle_pa) - np.sin(angle_pa) * y0
        y = y0 * np.cos(angle_pa) + np.sin(angle_pa) * x0


        ax.plot(x, y, label='Radius (arcsec)', color='black', linewidth=1)
        i +=1
def create_radius(angle_incl, angle_pa ):

    spacing = 10 #simple start -
    phi = np.linspace(0, 2*np.pi, spacing) # for now, 10
    #radius_plots = i  # arc sec for now #change this later for multiple radius

    x0 = 1 * np.cos(phi)
    y0 = 1 * np.cos(np.radians(angle_incl)) * np.sin(phi)

    # rotate - works for when east is to the left
    angle_pa = np.radians(angle_pa)
    x = x0 * np.cos(angle_pa) - np.sin(angle_pa) * y0
    y = y0 * np.cos(angle_pa) + np.sin(angle_pa) * x0
    return x, y

def plot_vect_on_radius(x_axis, y_axis,ax, bmaj, bmin , U_arr, Q_arr, pf):
    # sample points at 1/2 beams
    x_sample, y_sample = create_radius(34.9, 85.8)
    #plot_radius(34.9, 85.8, ax)


    # instantiate rectangular grid interpolator
    Q_interp = rgi((y_axis,x_axis), Q_arr, bounds_error=False, fill_value=np.nan)
    U_interp = rgi((y_axis,x_axis), U_arr, bounds_error=False, fill_value=np.nan)
    pf_interp = rgi((y_axis,x_axis), pf, bounds_error=False, fill_value=np.nan)

    points = np.column_stack((y_sample, x_sample))
    # interpolate from given points
    Q_sample = Q_interp(points)
    U_sample = U_interp(points)
    pf_sample = pf_interp(points)

    Q_sample = Q_sample.reshape(y_sample.shape)
    U_sample = U_sample.reshape(y_sample.shape)
    pf_sample = pf_sample.reshape(y_sample.shape)

    # calculate polarization angle
    pa_arr = .5 * np.arctan2(U_sample, Q_sample)  # this angle is east of north
    # cartesian
    # .3 arcsec for 1% polarization
    d_x = .3 * pf_sample * 100 * -np.sin(pa_arr)
    d_y = .3 * pf_sample * 100 * np.cos(pa_arr)

    # yy plotted opposite since ra is plotted on x-axis in conventional graph
    ax.quiver(y_sample,x_sample, d_y, d_x, pivot='middle', color='red', scale=1, scale_units='xy', headwidth=1e-10,
              headlength=1e-10, headaxislength=1e-10, width=0.005)