import numpy as np
import matplotlib.pyplot as plt
from radio_beam import Beam
import astropy.units as u



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

def plot_i(I, ra, dec, noise_I, **kwargs):
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
         artist = plot_beam(kwargs['bmaj'], kwargs['bmin'], kwargs['bpa'])
         ax.add_artist(artist)

    if 'ux' in kwargs:
        ax.quiver(ra, dec, kwargs['ux'], kwargs['uy'], pivot='middle', scale=30 )

    plt.xlim(3,-3) #manually swapped direction
    plt.ylim(-3,3)
    plt.show()


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

def plot_beam(bmaj, bmin,bpa):
    my_beam = Beam( bmaj * u.arcsec, bmin*u.arcsec, bpa*u.degree)
    y_cen_pix, x_cen_pix = -2,2
    pixscale = u.arcsec
    ellipse_artist = my_beam.ellipse_to_plot(x_cen_pix, y_cen_pix, pixscale, facecolor='black')
    return ellipse_artist
